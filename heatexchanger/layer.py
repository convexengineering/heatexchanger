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
    Q                [W]       heat transferred from  hot to cold fluid
    D_cold           [N]       total air drag
    D_hot            [N]       total water drag
    V_tot            [cm^3]    total volume
    V_mtrl           [cm^3]    volume of material
    g           9.81 [m*s^-2]  gravitational acceleration
    x_dim            [cm]      max hot length
    y_dim            [cm]      max cold length
    z_dim            [cm]      max height
    maxAR          5 [-]       max tile width variation
    T_max_hot        [K]       max temp. out
    T_min_cold       [K]       min temp. out
    T_in_hot         [K]       inlet temperature of hot fluid
    v_in_hot         [m/s]     inlet speed of hot fluid
    T_in_cold        [K]       inlet temperature of cold fluid
    v_in_cold        [m/s]     inlet speed of cold fluid
    solidity         [-]       solidity of HX
    max_solidity     [-]       max solidity allowed

    Upper Unbounded
    ---------------
    
    D_cold, D_hot, max_solidity, x_dim, y_dim, z_dim, T_max_hot, T_in_hot


    Lower Unbounded
    ---------------
    Q, T_in_cold

    """

    def setup(self, Ncoldpipes, Nhotpipes, coldfluid_model, hotfluid_model, material_model):
        self.Ncoldpipes = Ncoldpipes
        self.Nhotpipes = Nhotpipes

        exec parse_variables(Layer.__doc__)

        self.material = material_model()
        coldfluid = coldfluid_model()
        with Vectorize(Ncoldpipes):
            coldpipes = RectangularPipe(Nhotpipes, coldfluid,
                                        increasingT=True)
        self.coldpipes = coldpipes
        hotfluid = hotfluid_model()
        with Vectorize(Nhotpipes):
            hotpipes = RectangularPipe(Ncoldpipes, hotfluid,
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
            ("Hot_Channels", self.Nhotpipes),
            ("Cold_Channels", self.Ncoldpipes),
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

        # Unbounded variables from hierarchy
        self.P_in_hot = self.hotpipes.P_in
        self.P_in_cold = self.coldpipes.P_in
        self.P_out_hot = self.hotpipes.P_out
        self.P_out_cold = self.coldpipes.P_out
        self.T_in_hot = self.hotpipes.T_in
        self.T_in_cold = self.coldpipes.T_in

        with Vectorize(Nhotpipes):
            with Vectorize(Ncoldpipes):
                c = self.c = HXArea(self.material)

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
            c.x_cell >= c.n_fins_hot*(c.t_hot + hotpipes.w_fluid),
            c.y_cell >= c.n_fins_cld*(c.t_cld + coldpipes.w_fluid.T),
            # Making sure there is at least 1 (non-integer) number of fins
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
            T_min_cold <= T_max_hot,
            T_min_cold >= 1*units('K')
        ]

        for j in range(Nhotpipes):
            for i in range(Ncoldpipes):
                geom.extend([c.x_cell[i,j] == hotpipes.w[j],
                             c.y_cell[i,j] == coldpipes.w[i],
                             c.n_fins_hot[i,j] == hotpipes.n_fins[j],
                             c.n_fins_cld[i,j] == coldpipes.n_fins[i],
                             ])

        return [
            # SIZING
            geom,
            pipes,
            self.material,
            solidity == V_mtrl/V_tot,
            solidity <= max_solidity,

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
            V_mtrl >= (c.n_fins_hot*c.z_hot*c.t_hot*c.x_cell).sum()+(c.n_fins_cld*c.z_cld*c.t_cld*c.y_cell).sum()+(c.x_cell*c.y_cell*c.t_plate).sum(),
        ]
