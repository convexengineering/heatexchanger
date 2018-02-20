from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled, units
from gpkit.constraints.bounded import Bounded
from fluids import Air, Water
from hxarea import HXArea
from rectpipe import RectangularPipe
import numpy as np


class Layer(Model):
    """
    Variables
    ---------
    Q          [W]       heat transferred from air to liquid
    w          [m]       width
    d          [m]       depth
    h          [m]       height
    V        1 [cm^3]    volume
    h_liq      [m]       height of liquid layer
    h_air      [m]       height of air layer

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
                                       substitutions={"T_in": 303*units('K'),
                                                      "v_in": 20*units('m/s')})
            self.airpipes = airpipes
        water = Water()
        with Vectorize(Nwaterpipes):
            waterpipes = RectangularPipe(Nairpipes, water, increasingT=False,
                                         substitutions={"T_in": 500*units('K'),
                                                        "v_in": 5*units('m/s')})
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
                for j in range(Nairpipes):
                    waterCf.extend([waterpipes.l[i,j] >= sum(airpipes.w[0:j+1])])
            for i in range(Nairpipes):
                for j in range(Nwaterpipes):
                    airCf.extend([airpipes.l[i,j] >= sum(waterpipes.w[0:j+1])])

        return [
            # SIZING
            V >= d*w*h,
            h >= h_liq + h_air,
            waterpipes,
            d >= waterpipes.w.sum(),
            h_liq == waterpipes.h,
            airpipes,
            w >= airpipes.w.sum(),
            h_air == airpipes.h,

            # CONSERVATION OF HEAT
            SP_Qsum,
            c.dQ == waterpipes.dQ,
            c.dQ == airpipes.dQ.T,   # airpipes are rotated 90deg

            # HEAT EXCHANGE
            c, c.A_hx == airpipes.w.outer(waterpipes.w),
            c.T_hot == waterpipes.T[1:],  # Tcell = Tout (conservative)
            # NOTE: cell temperature could instead be a geometric mean
            #       of input and output
            c.T_cld == airpipes.T[1:].T,  # airpipes are rotated 90deg

            #DRAG
            #waterCf,
            #airCf,
        ]


if __name__ == "__main__":
    m = Layer(5, 5)
    m.cost = 1/m.Q + m.waterpipes.l.sum()*units('1/m/W') + m.airpipes.l.sum()*units('1/m/W') + m.waterpipes.Cf.sum()*units('1/W') + m.waterpipes.Cf.sum()*units('1/W')
    m = Model(m.cost, Bounded(m))
    sol = m.localsolve(verbosity=4)
    print sol(m.Q)
