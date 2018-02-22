from gpkit import Model, parse_variables


class RectangularPipe(Model):
    """
    Variables
    ---------
    mdot     [kg/s]   mass flow rate
    w        [m]      pipe width
    h        [m]      pipe height
    A        [m^2]    pipe area
    T_in     [K]      input temperature
    v_in     [m/s]    input velocity
    D        [N]      total drag

    Variables of length Nsegments+1
    -------------------------------
    v        [m/s]    fluid velocity
    T        [K]      fluid temperature

    Variables of length Nsegments
    -----------------------------
    dT       [K]      Change in fluid temperature over segment
    dQ       [W]      Magnitude of heat transfer over segment
    v_avg    [m/s]    Average fluid velocity over segment
    l        [m]      Reference flow length
    Re       [-]      Reynolds number
    Cf       [-]      Coefficient of friction over segment 
    D_seg    [N]      Drag over segment

    Upper Unbounded
    ---------------
    w, h, mdot, T_out (if increasingT), l

    Lower Unbounded
    ---------------
    dQ, T_out (if not increasingT)

    """
    def setup(self, Nsegments, fluid, increasingT):
        exec parse_variables(RectangularPipe.__doc__)
        self.increasingT = increasingT
        self.T_out = T[-1]
        if increasingT:
            temp = T[1:] >= T[:-1] + dT
        else:
            temp = T[:-1] >= T[1:] + dT
        return [
            fluid, temp,
            T[0] == T_in, v[0] == v_in,
            mdot == fluid.rho*v*A,
            A == w*h,
            dQ <= mdot*fluid.c*dT,
            v_avg**2 == v[0:-1]*v[1:],
            Re == (fluid.rho*v_avg*l/fluid.mu),
            Cf**5*Re == (0.059)**5,
        ]
