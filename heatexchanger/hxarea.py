from gpkit import Model, parse_variables, SignomialsEnabled

# Assumes hot water going S->N, and cold air going W->E, [0,0] at RH corner
class HXArea(Model):
    """
    Variables
    ---------
    h_hot       [W/K/m^2]    convective heat transfer coefficient of liquid
    h_cld      [W/K/m^2]    convective heat transfer coefficient of air
    T_hot       [K]          liquid temperature in cell
    T_cld       [K]          air temperature in cell
    T_r         [K]          wall temperature in cell
    A_hx        [m^2]        area of contact between liquid and air
    dQ          [W]          heat transferred from liquid to air by cell
    x_cell      [m]          x-width of cell
    y_cell      [m]          y-width of cell

    Upper Unbounded
    ---------------
    T_hot, T_r

    Lower Unbounded
    ---------------
    dQ, T_cld, T_r

    """
    def setup(self):
        exec parse_variables(HXArea.__doc__)
        with SignomialsEnabled():  # note that this turns into a posynomial
            return [dQ <= (T_hot-T_r)*h_hot*A_hx,
                    dQ <= (T_r-T_cld)*h_cld*A_hx,
                    A_hx == x_cell * y_cell]

