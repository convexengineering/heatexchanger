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
    D_wat      0.01 [N]       total water drag
    V_tot           [cm^3]    total volume
    V_mtrl          [cm^3]    volume of material
    g          9.81 [m*s^-2]  gravitational acceleration
    x_dim      5    [cm]      max hot length
    y_dim      10   [cm]      max cold length
    z_dim      1    [cm]      max height
    n_fins          [-]       fins per tile
    maxAR      5    [-]       max tile width variation
    T_max_hot  450  [K]       max temp. out
    T_min_cld    1  [K]       min temp. out
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
            airpipes = RectangularPipe(Nwaterpipes, n_fins, air,
                                       increasingT=True)
            self.airpipes = airpipes
        water = self.water = Water()
        with Vectorize(Nwaterpipes):
            waterpipes = RectangularPipe(Nairpipes, n_fins, water,
                                         increasingT=False)
            self.waterpipes = waterpipes
        pipes = [airpipes,
                 airpipes.T_in == T_in_air, airpipes.v_in == v_in_air,
                 waterpipes,
                 waterpipes.T_in == T_in_water, waterpipes.v_in == v_in_water]

        self.design_parameters = OrderedDict([
            # add T_max_hot, T_min_cld?
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
                c = self.c = HXArea(n_fins, self.material)

        waterCf = []
        airCf = []

        with SignomialsEnabled():
            # NOTE: unfortunately this appears unavoidable.
            #       perhaps an entropy-based approach could get around it?
            #       as the mass flows in each pipe become quite similar,
            #       it's already alllllmost GP, solving in 3-9 GP solves
            SP_Qsum = Q <= c.dQ.sum()

        geom = [
            V_tot >= waterpipes.V_seg.sum() + airpipes.V_seg.sum() + V_mtrl,
            c.x_cell == airpipes.l_seg.T,
            c.y_cell == waterpipes.l_seg,
            maxAR >= c.y_cell/c.x_cell,
            maxAR >= c.x_cell/c.y_cell,
            # Differentiating between flow width and cell width
            c.x_cell >= n_fins*(c.t_hot + waterpipes.w_fluid),
            c.y_cell >= n_fins*(c.t_cld + airpipes.w_fluid.T),
            # Making sure there is at least 1 (non-integer) number of fins
            n_fins >= 1.,
            c.dQ == waterpipes.dQ,
            c.dQ == airpipes.dQ.T,
            c.Tr_hot == waterpipes.Tr_int,
            c.Tr_cld == airpipes.Tr_int.T,
            c.T_hot == waterpipes.T_avg,
            c.T_cld == airpipes.T_avg.T,
            c.h_hot == waterpipes.h,
            c.h_cld == airpipes.h.T,
            c.z_hot == waterpipes.h_seg,
            c.z_cld == airpipes.h_seg.T,
            x_dim >= waterpipes.w.sum(),
            y_dim >= airpipes.w.sum(),
            z_dim >= c.z_hot+c.z_cld+c.t_plate,
            T_max_hot >= c.T_hot[-1, :],
            T_min_cld <= c.T_cld[:, -1],
            T_min_cld <= T_max_hot
        ]

        for j in range(Nwaterpipes):
            for i in range(Nairpipes):
                geom.extend([c.x_cell[i,j] == waterpipes.w[j],
                             c.y_cell[i,j] == airpipes.w[i],
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

            # TOTAL VOLUME REQUIREMENT
            V_tot <= x_dim*y_dim*z_dim,

            # MATERIAL VOLUME
            V_mtrl >= (n_fins*c.z_hot*c.t_hot*c.x_cell).sum()+(n_fins*c.z_cld*c.t_cld*c.y_cell).sum()+(c.x_cell*c.y_cell*c.t_plate).sum(),
        ]
