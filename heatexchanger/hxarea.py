from gpkit import Model, parse_variables, SignomialsEnabled, units

# Assumes hot water going S->N, and cold air going W->E, [0,0] at RH corner
class HXArea(Model):
    """
    Variables
    ---------
    dQ          [W]          heat transferred in cell
    h_hot       [W/K/m^2]    convective heat transfer coefficient of hot fluid
    h_cld       [W/K/m^2]    convective heat transfer coefficient of cold fluid
    T_hot       [K]          hot fluid temperature in cell
    T_cld       [K]          cold fluid temperature in cell
    T_r         [K]          wall internal temperature in cell
    Tr_hot      [K]          hot fluid-wall interface temperature
    Tr_cld      [K]          cold fluid-wall interface temperature
    x_cell      [cm]         x-width of cell
    y_cell      [cm]         y-width of cell
    z_hot       [cm]         z-height of hot cell
    z_cld       [cm]         z_height of cold cell
    t_plate     [cm]         thickness of separating plate
    t_hot       [cm]         fin thickness in hot fluid
    t_cld       [cm]         fin thickness in cold fluid
    A_hot       [cm^2]       area of hot fin
    A_cld       [cm^2]       area of cold fin
    n_fins_hot  [-]          number of fins on hot side
    n_fins_cld  [-]          number of fins on cold side


     Upper Unbounded
     ---------------
     T_hot, x_cell, y_cell, t_hot, t_cld, n_fins_hot, n_fins_cld

     Lower Unbounded
     ---------------
     dQ, T_cld

    """
    def setup(self, material):
        exec parse_variables(HXArea.__doc__)
        with SignomialsEnabled():  # note that these turn into posynomials
            dQ_definition = [dQ <= (T_hot-Tr_hot)*h_hot*A_hot,
                             dQ <= (Tr_cld-T_cld)*h_cld*A_cld]
        return [material, dQ_definition,
                A_hot == n_fins_hot*(2*y_cell*z_hot),
                A_cld == n_fins_cld*(2*x_cell*z_cld),
                Tr_hot   >= T_r + 0.33*(dQ/n_fins_hot*z_hot/(material.k*t_hot*y_cell)), #TODO: Refine
                T_r      >= Tr_cld + 0.33*(dQ/n_fins_cld*z_cld/(material.k*t_cld*x_cell)), #TODO: Refine
                t_plate  == material.t_min,
                t_hot    >= material.t_min,
                t_cld    >= material.t_min,
                ]
