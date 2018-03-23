from gpkit import Model, parse_variables, SignomialsEnabled, units

# Assumes hot water going S->N, and cold air going W->E, [0,0] at RH corner
class HXArea(Model):
    """
    Variables
    ---------
    dQ          [W]          heat transferred from liquid to air by cell
    h_hot       [W/K/m^2]    convective heat transfer coefficient of liquid
    h_cld       [W/K/m^2]    convective heat transfer coefficient of air
    T_hot       [K]          liquid temperature in cell
    T_cld       [K]          air temperature in cell
    T_r         [K]          wall internal temperature in cell
    Tr_hot      [K]          hot fluid-wall interface temperature
    Tr_cld      [K]          cold fluid-wall interface temperature
    x_cell      [m]          x-width of cell
    y_cell      [m]          y-width of cell
    z_hot       [m]          z-height of hot cell
    z_cld       [m]          z_height of cold cell
    t_plate     [m]          thickness of separating plate
    t_hot       [m]          fin thickness in liquid
    t_cld       [m]          fin thickness in gas
    A_hot       [m^2]        area of hot fin
    A_cld       [m^2]        area of cold fin

    """



    def setup(self, material):
        exec parse_variables(HXArea.__doc__)
        with SignomialsEnabled():  # note that this turns into a posynomial
            return [material,
                    dQ        <= (T_hot-Tr_hot)*h_hot*A_hot,
                    dQ        <= (Tr_cld-T_cld)*h_cld*A_cld,
                    A_hot**2  == 4*(x_cell * y_cell)*(2*y_cell*z_hot),
                    A_cld**2  == 4*(x_cell * y_cell)*(2*x_cell*z_cld),
                    T_hot     >= Tr_hot,
                    Tr_hot    >= T_r + 0.5*dQ*t_plate/(material.k*A_hot) + 0.001*(dQ*z_hot/(material.k*t_hot*y_cell)), #TODO: Refine
                    T_r       >= Tr_cld + 0.5*dQ*t_plate/(material.k*A_cld) + 0.001*(dQ*z_cld/(material.k*t_cld*x_cell)), #TODO: Refine
                    Tr_cld    >= T_cld,
                    t_plate   >= 0.01*units('cm'),
                    t_hot     >= 0.01*units('cm'),
                    t_cld     >= 0.01*units('cm'),
                    ]
