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
    Q               [W]       heat power from hot to cold fluid
    D_cold           [N]      total cold fluid drag
    D_hot           [N]       total hot fluid drag
    V_tot           [cm^3]    total volume
    V_mtrl          [cm^3]    volume of material
    g          9.81 [m*s^-2]  gravitational acceleration
    x_dim      5    [cm]      max hot length
    y_dim      10   [cm]      max cold length
    z_dim      1    [cm]      max height
    n_fins          [-]       fins per tile
    maxAR      5    [-]       max tile width variation
    T_max_hot  450  [K]       max temp. out
    T_max_cld       [K]       min temp. out

    Upper Unbounded
    ---------------
    D_cold, D_hot

    Lower Unbounded
    ---------------
    Q

    """
    def setup(self, Ncoldpipes, Nhotpipes, coldFluid, hotFluid, material):
        self.Ncoldpipes = Ncoldpipes
        self.Nhotpipes = Nhotpipes
        self.coldFluid = coldFluid
        self.hotFluid = hotFluid
        self.material = material

        exec parse_variables(Layer.__doc__)


        with Vectorize(Ncoldpipes):
            coldpipes = RectangularPipe(Nhotpipes, n_fins, coldFluid, increasingT=True,
                                       substitutions={"T_in": 303,
                                                      "v_in": 20})
            self.coldpipes = coldpipes
        with Vectorize(Nhotpipes):
            hotpipes = RectangularPipe(Ncoldpipes, n_fins, hotFluid, increasingT=False,
                                         substitutions={"T_in": 500,
                                                        "v_in": 1})
            self.hotpipes = hotpipes

        with Vectorize(Nhotpipes):
            with Vectorize(Ncoldpipes):
                c = self.c = HXArea(n_fins, self.material)

        hotCf = []
        coldCf = []

        with SignomialsEnabled():
            # NOTE: unfortunately this appears unavoidable.
            #       perhaps an entropy-based approach could get around it?
            #       as the mass flows in each pipe become quite similar,
            #       it's already alllllmost GP, solving in 3-9 GP solves
            SP_Qsum = Q <= c.dQ.sum()

        geom = [V_tot >= sum(sum(hotpipes.V_seg)) + sum(sum(coldpipes.V_seg)) + V_mtrl]

        # REMEMBER: LEFTMOST, INNERMOST RULE OF VECTORIZATION
        for i in range(Nhotpipes):
            for j in range(Ncoldpipes):
                print c.x_cell[j,i]
                geom.extend([c.x_cell[j,i] == hotpipes.w[i],
                             c.x_cell[j,i] == coldpipes.l_seg[i,j],
                             c.y_cell[j,i] == coldpipes.w[j],
                             c.y_cell[j,i] == hotpipes.l_seg[j,i],
                             maxAR >= c.y_cell[j,i]/c.x_cell[j,i],
                             maxAR >= c.x_cell[j,i]/c.y_cell[j,i],
                             # Arbitrary bounding for convergence.
                             c.x_cell[j,i] >= 0.1*units('cm'),
                             c.y_cell[j,i] >= 0.1*units('cm'),
                             # Differentiating between flow width and cell width
                             c.x_cell[j,i] >= n_fins*(c.t_hot[j,i] + hotpipes.w_fluid[j,i]),
                             c.y_cell[j,i] >= n_fins*(c.t_cld[j,i] + coldpipes.w_fluid[i,j]),
                             # Making sure there is at least 1 (non-integer) number of fins
                             n_fins >= 1.,
                             ])

        # Linking pipes in HXArea
        for i in range(Nhotpipes):
            for j in range(Ncoldpipes):
                geom.extend([
                        c.dQ[j,i]     == hotpipes.dQ[j,i],
                        c.dQ[j,i]     == coldpipes.dQ[i,j],
                        c.Tr_hot[j,i] == hotpipes.Tr_int[j,i],
                        c.Tr_cld[j,i] == coldpipes.Tr_int[i,j],
                        c.T_hot[j,i]  == hotpipes.T_avg[j,i],
                        c.T_cld[j,i]  == coldpipes.T_avg[i,j],
                        c.h_hot[j,i]  == hotpipes.h[j,i],
                        c.h_cld[j,i]  == coldpipes.h[i,j],
                        c.z_hot[j,i]  == hotpipes.h_seg[j,i],
                        c.z_cld[j,i]  == coldpipes.h_seg[i,j],
                        hotpipes.h_seg[j,i] >= 0.1*units('cm'),
                        hotpipes.h_seg[j,i] <= 1.0*units('cm'),
                        coldpipes.h_seg[i,j] >= 0.1*units('cm'),
                        coldpipes.h_seg[i,j] <= 1.0*units('cm'),
                        x_dim >= hotpipes.w.sum(),
                        y_dim >= coldpipes.w.sum(),
                        z_dim >= c.z_hot+c.z_cld+c.t_plate,
                        c.T_hot[-1,:] <= T_max_hot,
                        c.T_cld[:,-1] >= T_max_cld,
                        T_max_cld <= T_max_hot,
                        T_max_cld >= 1*units('K'),
                        ])

        return [
            # SIZING
            geom,
            hotpipes,
            coldpipes,
            self.material,

            # CONSERVATION OF HEAT
            SP_Qsum,
            c,

            #DRAG
            D_hot >= self.hotpipes.D.sum(),
            D_cold >= self.coldpipes.D.sum(),
            hotCf,
            coldCf,

            # TOTAL VOLUME REQUIREMENT
            V_tot <= 100*units('cm^3'),

            # MATERIAL VOLUME
            V_mtrl >= (n_fins*c.z_hot*c.t_hot*c.x_cell).sum()+(n_fins*c.z_cld*c.t_cld*c.y_cell).sum()+(c.x_cell*c.y_cell*c.t_plate).sum(),
        ]
