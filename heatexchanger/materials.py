from gpkit import Model, parse_variables


class Air(Model):
    """
    Variables
    ---------
    rho     1.2       [kg/m^3]    density of air
    c       1000      [J/K/kg]    heat capacity of air
    mu      1.81e-5   [Pa*s]      dynamic viscosity of air
    k       0.0262    [W/(m*K)]   thermal conductivity of air
    """
    def setup(self):
        exec parse_variables(Air.__doc__)


class Water(Model):
    """m.debug()
    Variables
    ---------
    rho     1000       [kg/m^3]    density of water
    c       4184       [J/K/kg]    heat capacity of water
    mu      8.90e-4    [Pa*s]      dynamic viscosity of water
    k       0.606      [W/(m*K)]   thermal conductivity of water
    """
    def setup(self):
        exec parse_variables(Water.__doc__)

class StainlessSteel(Model):
    """
    Variables
    ---------
    rho     7700       [kg/m^3]    density of steel
    c       500        [J/K/kg]    heat capacity of steel
    k       19         [W/(m*K)]   thermal conductivity of steel
    """
    def setup(self):
        exec parse_variables(Water.__doc__)
