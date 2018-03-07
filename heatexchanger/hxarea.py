from gpkit import Model, parse_variables, SignomialsEnabled


class HXArea(Model):
    """
    Variables
    ---------
    h_hx     10 [W/K/m^2]    coefficient of convection heat exchange
    T_hot       [K]          liquid temperature in cell
    T_cld       [K]          air temperature in cell
    T_r         [K]          wall temperature in cell
    A_hx        [m^2]        area of contact between liquid and air
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
            return [dQ <= (T_hot-T_cld)*h_hx*A_hx]
