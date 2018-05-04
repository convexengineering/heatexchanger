from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled
from materials import Air, Water, StainlessSteel
from hxarea import HXArea
from rectpipe import RectangularPipe
from collections import OrderedDict


class Layer(Model):
    """Combines heat exchanger pipes into a 2D layer

    Variables
    ---------
    Q                [W]       heat transferred from  hot to cold fluid
    D_cold      0.01 [N]       total air drag
    D_hot       0.01 [N]       total water drag
    V_tot            [cm^3]    total volume
    V_mtrl           [cm^3]    volume of material
    g           9.81 [m*s^-2]  gravitational acceleration
    x_dim          5 [cm]      max hot length
    y_dim         10 [cm]      max cold length
    z_dim          1 [cm]      max height
    n_fins           [-]       fins per tile
    maxAR          5 [-]       max tile width variation
    T_max_hot    450 [K]       max temp. out
    T_min_cold     1 [K]       min temp. out
    T_in_hot     500 [K]       inlet temperature of hot fluid
    v_in_hot       1 [m/s]     inlet speed of hot fluid
    T_in_cold    303 [K]       inlet temperature of cold fluid
    v_in_cold     20 [m/s]     inlet speed of cold fluid
    solidity         [-]       solidity of HX
    max_solidity 0.8 [-]       max solidity allowed

    Lower Unbounded
    ---------------
    Q

    """

    material_model = StainlessSteel
    coldfluid_model = Air
    hotfluid_model = Water

    def setup(self, Ncoldpipes, Nhotpipes):
        self.Ncoldpipes = Ncoldpipes
        self.Nhotpipes = Nhotpipes
        exec parse_variables(Layer.__doc__)

        self.material = self.material_model()
        with Vectorize(Nhotpipes):
            with Vectorize(Ncoldpipes):
                cells = self.cells = HXArea(n_fins, self.material)

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
        pipes = [
            coldpipes,
            coldpipes.T_in == T_in_cold, coldpipes.v_in == v_in_cold,
            hotpipes,
            hotpipes.T_in == T_in_hot, hotpipes.v_in == v_in_hot
        ]

        self.design_parameters = OrderedDict([
            # TODO: add T_max_hot, T_min_cold?
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
            ("max_solidity", max_solidity)
        ])

        geom = [
            V_tot >= hotpipes.V_seg.sum() + coldpipes.V_seg.sum() + V_mtrl,
            cells.x_cell == coldpipes.l_seg.T,
            cells.y_cell == hotpipes.l_seg,
            maxAR >= cells.y_cell/cells.x_cell,
            maxAR >= cells.x_cell/cells.y_cell,
            # Differentiating between flow width and cell width
            cells.x_cell >= n_fins*(cells.t_hot + hotpipes.w_fluid),
            cells.y_cell >= n_fins*(cells.t_cld + coldpipes.w_fluid.T),
            n_fins >= 1.,  # Making sure there is at least 1 fin
            cells.dQ == hotpipes.dQ,
            cells.dQ == coldpipes.dQ.T,
            cells.Tr_hot == hotpipes.Tr_int,
            cells.Tr_cld == coldpipes.Tr_int.T,
            cells.T_hot == hotpipes.T_avg,
            cells.T_cld == coldpipes.T_avg.T,
            cells.h_hot == hotpipes.h,
            cells.h_cld == coldpipes.h.T,
            cells.z_hot == hotpipes.h_seg,
            cells.z_cld == coldpipes.h_seg.T,
            x_dim >= hotpipes.w.sum(),
            y_dim >= coldpipes.w.sum(),
            z_dim >= cells.z_hot + cells.z_cld + cells.t_plate,
            T_max_hot >= cells.T_hot[-1, :],
            T_min_cold <= cells.T_cld[:, -1],
            T_min_cold <= T_max_hot
        ]

        for j in range(Nhotpipes):
            for i in range(Ncoldpipes):
                geom.extend([
                    cells.x_cell[i, j] == hotpipes.w[j],
                    cells.y_cell[i, j] == coldpipes.w[i],
                ])

        with SignomialsEnabled():
            SP_Qsum = Q <= cells.dQ.sum()

        return [
            SP_Qsum, cells, pipes, geom, self.material,

            # SOLIDITY
            solidity == V_mtrl/V_tot,
            solidity <= max_solidity,

            # DRAG
            D_hot >= self.hotpipes.D.sum(),
            D_cold >= self.coldpipes.D.sum(),

            # TOTAL VOLUME REQUIREMENT
            V_tot <= x_dim*y_dim*z_dim,

            # MATERIAL VOLUME
            V_mtrl >= ((n_fins*cells.z_hot*cells.t_hot*cells.x_cell).sum()
                       + (n_fins*cells.z_cld*cells.t_cld*cells.y_cell).sum()
                       + (cells.x_cell*cells.y_cell*cells.t_plate).sum()),
        ]
