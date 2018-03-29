from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled, units
from gpkit.constraints.bounded import Bounded
from materials import Air, Water, StainlessSteel
from hxarea import HXArea
from rectpipe import RectangularPipe
from relaxed_constants import relaxed_constants, post_process
import numpy as np


class Layer(Model):
    """
    Combines heat exchanger pipes into a 2D layer

    Variables
    ---------
    Q               [W]       heat transferred from air to liquid
    D_air           [N]       total air drag
    D_wat           [N]       total water drag
    V_tot           [cm^3]    total volume
    V_mtrl          [cm^3]    volume of material
    g          9.81 [m*s^-2]  gravitational acceleration
    x_dim      5    [cm]      max hot length
    y_dim      10   [cm]      max cold length
    z_dim      1    [cm]      max height
    maxAR      5    [-]       max tile width variation
    T_max_hot  450  [K]       max temp. out
    T_max_cld       [K]       min temp. out

    Upper Unbounded
    ---------------
    D_air, D_wat

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
            airpipes = RectangularPipe(Nwaterpipes, air, increasingT=True,
                                       substitutions={"T_in": 303,
                                                      "v_in": 20})
            self.airpipes = airpipes
        water = self.water = Water()
        with Vectorize(Nwaterpipes):
            waterpipes = RectangularPipe(Nairpipes, water, increasingT=False,
                                         substitutions={"T_in": 500,
                                                        "v_in": 1})
            self.waterpipes = waterpipes

        with Vectorize(Nwaterpipes):
            with Vectorize(Nairpipes):
                c = self.c = HXArea(self.material)

        waterCf = []
        airCf = []

        with SignomialsEnabled():
            # NOTE: unfortunately this appears unavoidable.
            #       perhaps an entropy-based approach could get around it?
            #       as the mass flows in each pipe become quite similar,
            #       it's already alllllmost GP, solving in 3-9 GP solves
            SP_Qsum = Q <= c.dQ.sum()

        geom = [V_tot >= sum(sum(waterpipes.V_seg)) + sum(sum(airpipes.V_seg)) + V_mtrl]

        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                print c.x_cell[i,j]
                geom.extend([c.x_cell[i,j] == waterpipes.w[i],
                             c.x_cell[i,j] == airpipes.l_seg[j,i],
                             c.y_cell[i,j] == airpipes.w[j],
                             c.y_cell[i,j] == waterpipes.l_seg[i,j],
                             maxAR >= c.y_cell[i,j]/c.x_cell[i,j],
                             maxAR >= c.x_cell[i,j]/c.y_cell[i,j],
                             # Arbitrary bounding for convergence.
                             c.x_cell[i,j] >= 0.1*units('cm'),
                             c.y_cell[i,j] >= 0.1*units('cm'),
                             ])

        # Linking pipes in c
        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                geom.extend([
                        c.dQ[i,j]     == waterpipes.dQ[i,j],
                        c.dQ[i,j]     == airpipes.dQ[j,i],
                        c.Tr_hot[i,j] == waterpipes.Tr_int[i,j],
                        c.Tr_cld[i,j] == airpipes.Tr_int[j,i],
                        c.T_hot[i,j]  == waterpipes.T_avg[i,j],
                        c.T_cld[i,j]  == airpipes.T_avg[j,i],
                        c.h_hot[i,j]  == waterpipes.h[i,j],
                        c.h_cld[i,j]  == airpipes.h[j,i],
                        c.z_hot[i,j]  == waterpipes.h_seg[i,j],
                        c.z_cld[i,j]  == airpipes.h_seg[j,i],
                        c.t_hot[i,j]     >= 0.05/((i+1.)**3*(j+1.)**3.)**(1./3.)*units('cm'),
                        c.t_cld[i,j]     >= 0.05/((i+1.)**3*(j+1.)**3.)**(1./3.)*units('cm'),
                        waterpipes.h_seg[i,j] >= 0.1*units('cm'),
                        waterpipes.h_seg[i,j] <= 1.0*units('cm'),
                        airpipes.h_seg[j,i] >= 0.1*units('cm'),
                        airpipes.h_seg[j,i] <= 1.0*units('cm'),
                        x_dim >= waterpipes.w.sum(),
                        y_dim >= airpipes.w.sum(),
                        z_dim >= c.z_hot+c.z_cld+c.t_plate,
                        c.T_hot[-1,:] <= T_max_hot,
                        c.T_cld[:,-1] >= T_max_cld,
                        T_max_cld <= T_max_hot,
                        T_max_cld >= 1*units('K'),
                        ])

        return [
            # SIZING
            geom,
            waterpipes,
            airpipes,
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
    m = Layer(5, 6)
    m.cost = (m.D_air+m.D_wat)/m.Q
    sol = m.localsolve(verbosity=2)
    print sol("Q")

    with open("sol.txt", "w") as f:
        f.write(sol.table())

    from writetotext import genHXData
    genHXData(m, sol)
