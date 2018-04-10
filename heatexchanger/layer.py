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
    D_cold     0.01 [N]       total air drag
    D_hot      0.01 [N]       total water drag
    V_tot           [cm^3]    total volume
    V_mtrl          [cm^3]    volume of material
    g          9.81 [m*s^-2]  gravitational acceleration
    x_dim      5    [cm]      max hot length
    y_dim      10   [cm]      max cold length
    z_dim      1    [cm]      max height
    n_fins          [-]       fins per tile
    maxAR      5    [-]       max tile width variation
    T_max_hot  450  [K]       max temp. out
    T_min_cold   1  [K]       min temp. out
    T_in_hot 500    [K]       inlet temperature of hot fluid
    v_in_hot 1      [m/s]     inlet speed of hot fluid
    T_in_cold   303 [K]       inlet temperature of cold fluid
    v_in_cold   20  [m/s]     inlet speed of cold fluid

    Lower Unbounded
    ---------------
    Q

    """

    coldfluid_model = Air
    hotfluid_model = Water
    material_model = StainlessSteel

    def setup(self, Ncoldpipes, Nhotpipes):
        self.Ncoldpipes = Ncoldpipes
        self.Nhotpipes = Nhotpipes

        exec parse_variables(Layer.__doc__)

        self.material = self.material_model()
        coldfluid = self.coldfluid_model()
        with Vectorize(Ncoldpipes):
            coldpipes = RectangularPipe(Nhotpipes, n_fins, coldfluid,
                                        increasingT=True)
            self.coldpipes = coldpipes
        hotfluid = self.hotfluid_model()
        with Vectorize(Nhotpipes):
            hotpipes = RectangularPipe(Ncoldpipes, n_fins, hotfluid,
                                       increasingT=False)
            self.hotpipes = hotpipes
        pipes = [coldpipes,
                 coldpipes.T_in == T_in_cold, coldpipes.v_in == v_in_cold,
                 hotpipes,
                 hotpipes.T_in == T_in_hot, hotpipes.v_in == v_in_hot]

        self.design_parameters = OrderedDict([
            # add T_max_hot, T_min_cold?
            ("gravity", g),
            ("x_width", x_dim),
            ("y_width", y_dim),
            ("z_width", z_dim),
            ("Hot_Channels", self.Ncoldpipes),
            ("Cold_Channels", self.Nhotpipes),
            ("Hot_Drag", D_hot),
            ("Cold_Drag", D_cold),
            ("c_metal", self.material.c),
            ("k_metal", self.material.k),
            ("rho_metal", self.material.rho),
            ("t_min_metal", self.material.t_min),
            ("c_coldfluid", coldfluid.c),
            ("k_coldfluid", coldfluid.k),
            ("rho_coldfluid", coldfluid.rho),
            ("mu_coldfluid", coldfluid.mu),
            ("Ti_coldfluid", T_in_cold),
            ("vi_coldfluid", v_in_cold),
            ("c_hotfluid", hotfluid.c),
            ("k_hotfluid", hotfluid.k),
            ("rho_hotfluid", hotfluid.rho),
            ("mu_hotfluid", hotfluid.mu),
            ("Ti_hotfluid", T_in_hot),
            ("vi_hotfluid", v_in_hot),
        ])

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

        geom = [
            V_tot >= hotpipes.V_seg.sum() + coldpipes.V_seg.sum() + V_mtrl,
            c.x_cell == coldpipes.l_seg.T,
            c.y_cell == hotpipes.l_seg,
            maxAR >= c.y_cell/c.x_cell,
            maxAR >= c.x_cell/c.y_cell,
            # Differentiating between flow width and cell width
            c.x_cell >= n_fins*(c.t_hot + hotpipes.w_fluid),
            c.y_cell >= n_fins*(c.t_cld + coldpipes.w_fluid.T),
            # Making sure there is at least 1 (non-integer) number of fins
            n_fins >= 1.,
            c.dQ == hotpipes.dQ,
            c.dQ == coldpipes.dQ.T,
            c.Tr_hot == hotpipes.Tr_int,
            c.Tr_cld == coldpipes.Tr_int.T,
            c.T_hot == hotpipes.T_avg,
            c.T_cld == coldpipes.T_avg.T,
            c.h_hot == hotpipes.h,
            c.h_cld == coldpipes.h.T,
            c.z_hot == hotpipes.h_seg,
            c.z_cld == coldpipes.h_seg.T,
            x_dim >= hotpipes.w.sum(),
            y_dim >= coldpipes.w.sum(),
            z_dim >= c.z_hot+c.z_cld+c.t_plate,
            T_max_hot >= c.T_hot[-1, :],
            T_min_cold <= c.T_cld[:, -1],
            T_min_cold <= T_max_hot
        ]

        for j in range(Nhotpipes):
            for i in range(Ncoldpipes):
                geom.extend([c.x_cell[i,j] == hotpipes.w[j],
                             c.y_cell[i,j] == coldpipes.w[i],
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
            D_hot >= self.hotpipes.D.sum(),
            D_cold >= self.coldpipes.D.sum(),
            hotCf,
            coldCf,

            # TOTAL VOLUME REQUIREMENT
            V_tot <= x_dim*y_dim*z_dim,

            # MATERIAL VOLUME
            V_mtrl >= (n_fins*c.z_hot*c.t_hot*c.x_cell).sum()+(n_fins*c.z_cld*c.t_cld*c.y_cell).sum()+(c.x_cell*c.y_cell*c.t_plate).sum(),
        ]
