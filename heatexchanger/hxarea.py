from gpkit import Model, parse_variables, SignomialsEnabled, units

# Assumes hot water going S->N, and cold air going W->E, [0,0] at RH corner
class HXArea(Model):
    """
    Variables
    ---------
    h_hot       [W/K/m^2]    convective heat transfer coefficient of liquid
    h_cld       [W/K/m^2]    convective heat transfer coefficient of air
    T_hot       [K]          liquid temperature in cell
    T_cld       [K]          air temperature in cell
    T_r         [K]          wall internal temperature in cell
    Tr_hot      [K]          hot fluid-wall interface temperature
    Tr_cld      [K]          cold fluid-wall interface temperature
    A_hx        [m^2]        area of contact between liquid and air
    dQ          [W]          heat transferred from liquid to air by cell
    x_cell      [m]          x-width of cell
    y_cell      [m]          y-width of cell
    t_plate     [m]          thickness of separating plate

    Upper Unbounded
    ---------------
    T_hot, T_r

    Lower Unbounded
    ---------------
    dQ, T_cld, T_r

    """

    # t_hot       [m]          fin thickness in liquid
    # t_cld       [m]          fin thickness in gas
    # mL_hot      [-]          hot fin parameter
    # mL_cld      [-]          cold fin parameter
    # eta_hot     [-]          hot fin efficiency
    # eta_cld     [-]          cold fin efficiency 

    def setup(self, material):
        exec parse_variables(HXArea.__doc__)
        with SignomialsEnabled():  # note that this turns into a posynomial
            return [dQ     <= (T_hot-Tr_hot)*h_hot*A_hx,
                    dQ     <= (Tr_cld-T_cld)*h_cld*A_hx,
                    A_hx   == x_cell * y_cell,
                    T_hot  >= Tr_hot,
                    Tr_hot >= T_r + 0.5*dQ*t_plate/(material.k*A_hx),
                    T_r    >= Tr_cld + 0.5*dQ*t_plate/(material.k*A_hx),
                    Tr_cld >= T_cld,
                    t_plate >= 0.01*units('cm'),
                    #mL_cld == h_cld*(h_cld*x_cell)/(material.k*)
                    ]

