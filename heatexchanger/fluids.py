from gpkit import Model, parse_variables


class Air(Model):
    """
    Variables
    ---------
    rho     1.2       [kg/m^3]    density of air
    c       1         [J/K/kg]    heat capacity of air
    mu      1.81e-5   [Pa*s]        dynamic viscosity of air
    """
    def setup(self):
        exec parse_variables(Air.__doc__)


class Water(Model):
    """
    Variables
    ---------
    rho     1e3       [kg/m^3]    density of water
    c       3         [J/K/kg]    heat capacity of water
    mu      8.90e-4   [Pa*s]        dynamic viscosity of water
    """
    def setup(self):
        exec parse_variables(Water.__doc__)
