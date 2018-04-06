from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled, units
from gpkit.constraints.bounded import Bounded
from materials import Air, Water, StainlessSteel
from hxarea import HXArea
from rectpipe import RectangularPipe
from relaxed_constants import relaxed_constants, post_process
import numpy as np
from collections import OrderedDict


class Layer(Model):
    """
    Combines heat exchanger pipes into a 2D layer

    Variables
    ---------
    Q               [W]       heat transferred from air to liquid
    D_air      0.01 [N]       total air drag
    D_wat      0.1  [N]       total water drag
    V_tot           [cm^3]    total volume
    V_mtrl          [cm^3]    volume of material
    g          9.81 [m*s^-2]  gravitational acceleration
    x_dim      5    [cm]      max hot length
    y_dim      10    [cm]      max cold length
    z_dim      1    [cm]      max height
    maxAR      8    [-]       max aspect ratio of tiles
    T_in_water 500  [K]       inlet temperature of water
    v_in_water 1    [m/s]     inlet speed of water
    T_in_air   303  [K]       inlet temperature of air
    v_in_air   20   [m/s]     inlet speed of air

    Lower Unbounded
    ---------------
    Q

    """
    def setup(self, Nairpipes, Nwaterpipes):
        self.Nwaterpipes = Nwaterpipes
        self.Nairpipes = Nairpipes
        exec parse_variables(Layer.__doc__)

        self.material = StainlessSteel()

        air = self.air = Air()
        with Vectorize(Nairpipes):
            airpipes = RectangularPipe(Nwaterpipes, air, increasingT=True)
            self.airpipes = airpipes
        water = self.water = Water()
        with Vectorize(Nwaterpipes):
            waterpipes = RectangularPipe(Nairpipes, water, increasingT=False)
            self.waterpipes = waterpipes
        pipes = [airpipes,
                 airpipes.T_in == T_in_air, airpipes.v_in == v_in_air,
                 waterpipes,
                 waterpipes.T_in == T_in_water, waterpipes.v_in == v_in_water]

        self.design_parameters = OrderedDict([
            ("gravity", g),
            ("x_width", x_dim),
            ("y_width", y_dim),
            ("z_width", z_dim),
            ("Air_Channels", self.Nairpipes),
            ("Water_Channels", self.Nwaterpipes),
            ("Air_Drag", D_air),
            ("Water_Drag", D_wat),
            ("c_metal", self.material.c),
            ("k_metal", self.material.k),
            ("rho_metal", self.material.rho),
            ("t_min_metal", self.material.t_min),
            ("c_air", air.c),
            ("k_air", air.k),
            ("rho_air", air.rho),
            ("mu_air", air.mu),
            ("Ti_air", T_in_air),
            ("vi_air", v_in_air),
            ("c_water", water.c),
            ("k_water", water.k),
            ("rho_water", water.rho),
            ("mu_water", water.mu),
            ("Ti_water", T_in_water),
            ("vi_water", v_in_water),
        ])

        with Vectorize(Nwaterpipes):
            with Vectorize(Nairpipes):
                c = self.c = HXArea(self.material)

        waterCf = [waterpipes.D[i] >= waterpipes.fr[i]*waterpipes.A_seg[i, 0]
                   for i in range(Nwaterpipes)]
        airCf = [airpipes.D[i] >= airpipes.fr[i]*airpipes.A_seg[i, 0]
                 for i in range(Nairpipes)]

        with SignomialsEnabled():
            SP_Qsum = [Q <= c.dQ.sum()]

        geom = [
            V_tot >= waterpipes.V_seg.sum() + airpipes.V_seg.sum() + V_mtrl,
            maxAR >= c.x_cell/c.y_cell,
            maxAR >= c.y_cell/c.x_cell,
            c.dQ     == waterpipes.dQ,
            c.dQ     == airpipes.dQ.T,
            c.Tr_hot == waterpipes.Tr_int,
            c.Tr_cld == airpipes.Tr_int.T,
            c.T_hot  == waterpipes.T_avg,
            c.T_cld  == airpipes.T_avg.T,
            c.h_hot  == waterpipes.h,
            c.h_cld  == airpipes.h.T,
            c.z_hot  == waterpipes.h_seg,
            c.z_cld  == airpipes.h_seg.T,
            x_dim >= waterpipes.w.sum(),
            y_dim >= airpipes.w.sum()
        ]

        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                waterCf.extend([waterpipes.l[i,j] == (j+1)*airpipes.w[0:j+1].prod()**(1./(j+1))])
                airCf.extend([airpipes.l[j,i] == (i+1)*waterpipes.w[0:i+1].prod()**(1./(i+1))])

                geom.extend([
                        c.x_cell[i,j] == waterpipes.w[i],
                        c.x_cell[i,j] == airpipes.l_seg[j,i],
                        c.y_cell[i,j] == airpipes.w[j],
                        c.y_cell[i,j] == waterpipes.l_seg[i,j],
                        c.t_hot[i,j]  >= 0.05/((i+1.)**3*(j+1.)**3.)**(1./3.)*units('cm'),
                        c.t_cld[i,j]  >= 0.05/((i+1.)**3*(j+1.)**3.)**(1./3.)*units('cm'),
                    ])

        return [
            # SIZING
            geom,
            pipes,
            self.material,

            # CONSERVATION OF HEAT
            SP_Qsum,
            c,

            #DRAG
            D_wat >= self.waterpipes.D.sum(),
            D_air >= self.airpipes.D.sum(),
            waterCf,
            airCf,

            # HEAT EXCHANGE REQUIREMENT
            #waterpipes.T[:,-1] <= 400*units('K'),

            # TOTAL VOLUME REQUIREMENT
            V_tot <= 100*units('cm^3'),

            # MATERIAL VOLUME
            V_mtrl >= (c.z_hot*c.t_hot*c.x_cell).sum()+(c.z_cld*c.t_cld*c.y_cell).sum()+(c.x_cell*c.y_cell*c.t_plate).sum(),
        ]


if __name__ == "__main__":
    m = Layer(5, 5)
    m.cost = (m.D_air+m.D_wat)/m.Q
    sol = m.localsolve(verbosity=1)
    print sol("Q")

    with open("sol.txt", "w") as f:
        f.write(sol.table())

    from writetotext import genHXData
    genHXData(m, sol)
