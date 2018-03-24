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


    Upper Unbounded
    ---------------
    D_air, D_wat

    Lower Unbounded
    ---------------
    Q, airpipes.dP_scale, waterpipes.dP_scale

    """
    def setup(self, Nairpipes, Nwaterpipes):
        exec parse_variables(Layer.__doc__)
        self.Nairpipes = Nairpipes
        self.Nwaterpipes = Nwaterpipes
        self.material = StainlessSteel()

        air = Air()
        with Vectorize(Nairpipes):
            airpipes = RectangularPipe(Nwaterpipes, air, increasingT=True,
                                       substitutions={"T_in": 303,
                                                      "v_in": 20})
            self.airpipes = airpipes
        water = Water()
        with Vectorize(Nwaterpipes):
            waterpipes = RectangularPipe(Nairpipes, water, increasingT=False,
                                         substitutions={"T_in": 500,
                                                        "v_in": 5})
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
            for i in range(Nwaterpipes):
                waterCf.extend([waterpipes.D[i] >= waterpipes.fr[i]*waterpipes.A_seg[i,0],
                                waterpipes.fr[i] >= waterpipes.dP[i,:].sum(),
                                #waterpipes.dP_scale[0] == waterpipes.dP_scale[i]
                                ])
                for j in range(Nairpipes):
                    waterCf.extend([waterpipes.l[i,j] <= sum(airpipes.w[0:j+1]),
                                    waterpipes.dP[i,j] >= waterpipes.dP_scale[i]*0.5*water.rho*waterpipes.v_avg[i,j]**2*waterpipes.Cf[i,j]*waterpipes.l_seg[i,j]/waterpipes.dh[i,j],
                                            ])
            for i in range(Nairpipes):
                airCf.extend([airpipes.D[i] >= airpipes.fr[i]*airpipes.A_seg[i,0],
                              airpipes.fr[i] >= airpipes.dP[i,:].sum(),
                              #airpipes.dP_scale[0] == airpipes.dP_scale[i]
                              ])
                for j in range(Nwaterpipes):
                    airCf.extend([airpipes.l[i,j] <= sum(waterpipes.w[0:j+1]),
                                  airpipes.dP[i,j] >= airpipes.dP_scale[i]*0.5*air.rho*airpipes.v_avg[i,j]**2*airpipes.Cf[i,j]*airpipes.l_seg[i,j]/airpipes.dh[i,j],
                                            ])

        geom = [V_tot >= sum(sum(waterpipes.V_seg)) + sum(sum(airpipes.V_seg)) + V_mtrl]

        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                geom.extend([c.x_cell[i,j] == waterpipes.w[i],
                             c.x_cell[i,j] == airpipes.l_seg[j,i],
                             c.y_cell[i,j] == airpipes.w[j],
                             c.y_cell[i,j] == waterpipes.l_seg[i,j],
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
                        c.t_plate[i,j]   >= 0.01*units('cm'),
                        c.t_hot[i,j]     >= 0.01*units('cm'),
                        c.t_cld[i,j]     >= 0.01*units('cm'),
                        c.t_hot[i,j]     >= 0.05/((i+1.)**3*(j+1.)**3.)**(1./3.)*units('cm'),
                        c.t_cld[i,j]     >= 0.05/((i+1.)**3*(j+1.)**3.)**(1./3.)*units('cm'),
                        waterpipes.h_seg[i,j] >= 0.1*units('cm'),
                        waterpipes.h_seg[i,j] <= 0.5*units('cm'),
                        airpipes.h_seg[j,i] >= 0.1*units('cm'),
                        airpipes.h_seg[j,i] <= 0.5*units('cm'),
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

            # TOTAL VOLUME REQUIREMENT
            V_tot <= 2*units('cm^3'),

            # MATERIAL VOLUME
            V_mtrl >= (c.z_hot*c.t_hot*c.x_cell).sum()+(c.z_cld*c.t_cld*c.y_cell).sum()+(c.x_cell*c.y_cell*c.t_plate).sum(),
        ]


if __name__ == "__main__":
    m = Layer(5,5)
    # m.substitutions.update({
    #     'V_tot':1*units('cm^3'),
    #     'Q'    :4*units('W')
    #     })
    penalties = (m.waterpipes.dP_scale.prod()*m.airpipes.dP_scale.prod()*m.waterpipes.dT.prod()*m.airpipes.dT.prod())**-1
    m.cost = penalties*(m.D_air+m.D_wat)/m.Q
    #m = Model(m.cost,Bounded(m))
    #m = relaxed_constants(m)
    sol = m.localsolve(verbosity=4)
    #post_process(sol)
    print sol('Q')
    # print sol("eta_h")
