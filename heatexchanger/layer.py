from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled, units
from gpkit.constraints.bounded import Bounded
from fluids import Air, Water
from hxarea import HXArea
from rectpipe import RectangularPipe
from relaxed_constants import relaxed_constants
import numpy as np


class Layer(Model):
    """
    Variables
    ---------
    Q          [W]       heat transferred from air to liquid
    V_tot    1 [cm^3]    total volume

    Lower Unbounded
    ---------------
    Q

    """
    def setup(self, Nairpipes, Nwaterpipes):
        exec parse_variables(Layer.__doc__)
        self.Nairpipes = Nairpipes
        self.Nwaterpipes = Nwaterpipes

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
                c = self.c = HXArea()

        waterCf = []
        airCf = []

        with SignomialsEnabled():
            # NOTE: unfortunately this appears unavoidable.
            #       perhaps an entropy-based approach could get around it?
            #       as the mass flows in each pipe become quite similar,
            #       it's already alllllmost GP, solving in 3-9 GP solves
            SP_Qsum = Q <= c.dQ.sum()
            for i in range(Nwaterpipes):
                waterCf.extend([waterpipes.D[i] >= waterpipes.fr[i]*waterpipes.A_seg[i,0]])
                for j in range(Nairpipes):
                    waterCf.extend([waterpipes.l[i,j] <= sum(airpipes.w[0:j+1]),
                                    waterpipes.dP[i,j] >= 0.5*water.rho*waterpipes.v_avg[i,j]**2*waterpipes.Cf[i,j]*waterpipes.l_seg[i,j]/waterpipes.dh[i,j],
                                    waterpipes.D_seg[i,j] == 0.5*water.rho*waterpipes.v_avg[i,j]**2*waterpipes.Cf[i,j]*waterpipes.w[i]*waterpipes.l_seg[i,j],
                                            ])
            for i in range(Nairpipes):
                airCf.extend([airpipes.D[i] >= airpipes.fr[i]*airpipes.A_seg[i,0]])
                for j in range(Nwaterpipes):
                    airCf.extend([airpipes.l[i,j] <= sum(waterpipes.w[0:j+1]),
                                  airpipes.dP[i,j] >= 0.5*air.rho*airpipes.v_avg[i,j]**2*airpipes.Cf[i,j]*airpipes.l_seg[i,j]/airpipes.dh[i,j],
                                  airpipes.D_seg[i,j] == 0.5*air.rho*airpipes.v_avg[i,j]**2*airpipes.Cf[i,j]*airpipes.w[i]*airpipes.l_seg[i,j],
                                            ])

        geom = [V_tot >= sum(sum(waterpipes.V_seg)) + sum(sum(airpipes.V_seg))]

        # Arbitrary bounding for convergence.
        for i in range(Nwaterpipes):
            continue
            geom.extend([waterpipes.w[i] >= 0.1*units('cm'),
                         waterpipes.h_seg[i,:] >= 0.1*units('cm'),
                         ])
        for j in range(Nairpipes):
            continue
            geom.extend([airpipes.w[j] >= 0.1*units('cm'),
                         airpipes.h_seg[j,:] >= 0.1*units('cm'),
                         ])
        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                geom.extend([c.x_cell[i,j] == waterpipes.w[i],
                             c.x_cell[i,j] == airpipes.l_seg[j,i],
                             c.y_cell[i,j] == airpipes.w[j],
                             c.y_cell[i,j] == waterpipes.l_seg[i,j],
                             ])

        # Linking pipes in c
        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                geom.extend([
                        c.dQ[i,j] == waterpipes.dQ[i,j],
                        c.dQ[i,j] == airpipes.dQ[j,i],
                        c.T_r[i,j] == waterpipes.Tr_int[i,j],
                        c.T_r[i,j] == airpipes.Tr_int[j,i],
                        c.T_hot[i,j] == waterpipes.T_avg[i,j],
                        c.T_cld[i,j] == airpipes.T_avg[j,i],
                        c.T_hot[i,j] >= c.T_r[i,j],
                        c.T_r[i,j] >= c.T_cld[i,j],
                    ])

        return [
            # SIZING
            geom,
            waterpipes,
            airpipes,

            # CONSERVATION OF HEAT
            SP_Qsum,
            c,

            #DRAG
            waterCf,
            airCf,
        ]


if __name__ == "__main__":
    m = Layer(5, 5)
    m.cost = 1/m.Q + 100*m.waterpipes.D.sum()*units('1/(N*W)')+ 1000*m.airpipes.D.sum()*units('1/(N*W)') + 1000*m.waterpipes.D_seg.sum()*units('1/(N*W)')+ 100*m.airpipes.D_seg.sum()*units('1/(N*W)')
    #m = Model(m.cost,Bounded(m))
    m = relaxed_constants(m)
    sol = m.localsolve(verbosity=4)
    print sol('Q')
