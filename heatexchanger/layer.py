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
                waterCf.extend([waterpipes.D[i] >= waterpipes.fr[i]*waterpipes.w[i]*waterpipes.h_seg[i,0]])
                for j in range(Nairpipes):
                    waterCf.extend([waterpipes.l[i,j] <= sum(airpipes.w[0:j+1]),
                                    waterpipes.dP[i,j] >= 0.5*water.rho*waterpipes.v_avg[i,j]**2*waterpipes.Cf[i,j]*airpipes.w[j]/waterpipes.dh[i],
                                    waterpipes.D_seg[i,j] == 0.5*water.rho*waterpipes.v_avg[i,j]**2*waterpipes.Cf[i,j]*waterpipes.w[i]*airpipes.w[j],
                                            ])
            for i in range(Nairpipes):
                airCf.extend([airpipes.D[i] >= airpipes.fr[i]*airpipes.w[i]*airpipes.h_seg[i,0]])
                for j in range(Nwaterpipes):
                    airCf.extend([airpipes.l[i,j] <= sum(waterpipes.w[0:j+1]),
                                airpipes.dP[i,j] == 0.5*air.rho*airpipes.v_avg[i,j]**2*airpipes.Cf[i,j]*waterpipes.w[j]/airpipes.dh[i],
                                airpipes.D_seg[i,j] == 0.5*air.rho*airpipes.v_avg[i,j]**2*airpipes.Cf[i,j]*airpipes.w[i]*waterpipes.w[j],
                                            ])
        geom = [V_tot >= sum(sum(c.V_hot)) + sum(sum(c.V_cold))]
        for i in range(Nwaterpipes):
            for j in range(Nairpipes):
                    geom.extend([c.w_cell[i,j] == waterpipes.w[i],
                                 c.d_cell[i,j] == airpipes.w[j],
                                 c.w_cell[i,j] >= 1*units('cm'),
                                 c.d_cell[i,j] >= 1*units('cm'),
                                 c.h_cold[i,j] >= 0.1*units('cm'),
                                 c.h_hot[i,j] >= 0.1*units('cm')])
        return [
            # SIZING
            geom,
            waterpipes,
            airpipes,

            # CONSERVATION OF HEAT
            SP_Qsum,
            c.dQ == waterpipes.dQ,
            c.dQ == airpipes.dQ.T,   # airpipes are rotated 90deg
            c.T_r == waterpipes.Tr_int,
            c.T_r == airpipes.Tr_int.T,

            # HEAT EXCHANGE
            c, c.A_hx == airpipes.w.outer(waterpipes.w),
            c.T_hot == waterpipes.T[1:],  # Tcell = Tout (conservative)
            # NOTE: cell temperature could instead be a geometric mean
            #       of input and output
            c.T_cld == airpipes.T[1:].T,  # airpipes are rotated 90deg

            # Coupling geometry
            c.V_cold == waterpipes.V_seg,
            c.V_hot == airpipes.V_seg.T,
            c.h_hot == waterpipes.h_seg,
            c.h_cold == airpipes.h_seg.T,

            #DRAG
            waterCf,
            airCf,
        ]


if __name__ == "__main__":
    m = Layer(5, 5)
    m.cost = 1/m.Q + 100*m.waterpipes.D.sum()*units('1/(N*W)')+ 100*m.airpipes.D.sum()*units('1/(N*W)')
    #m = Model(m.cost,Bounded(m))
    m = relaxed_constants(m)
    sol = m.localsolve(verbosity=4)
    print sol('Q')
