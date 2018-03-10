from gpkit import Model, parse_variables, SignomialsEnabled

# Assumes hot water going S->N, and cold air going W->E
class HXArea(Model):
    """
    Variables
    ---------
    h_hx     10 [W/K/m^2]    coefficient of convection heat exchange
    T_hot       [K]          liquid temperature in cell
    T_cld       [K]          air temperature in cell
    T_r         [K]          wall temperature in cell
    A_hx        [m^2]        area of contact between liquid and air
    V_cold      [m^3]        volume of cold cell
    V_hot       [m^3]        volume of hot cell
    w_cell      [m]          width of cell in air direction
    d_cell      [m]          width of cell in water direction 
    h_cold      [m]          height of cold cell
    h_hot       [m]          height of hot cell
    dQ          [W]          heat transferred from liquid to air by cell


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
            return [dQ <= (T_hot-T_cld)*h_hx*A_hx,
                    V_cold == h_cold*d_cell*w_cell,
                    V_hot  == h_hot*d_cell*w_cell]

