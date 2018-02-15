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

    Variables of length Nsegments+1
    -------------------------------
    v        [m/s]    fluid velocity
    T        [K]      fluid temperature

    Variables of length Nsegments
    -----------------------------
    dT       [K]      Change in fluid temperature over segment
    dQ       [W]      Magnitude of heat transfer over segment
    Cf       [-]      Coefficient of friction over segment 
    l        [m]      Reference flow length

    Upper Unbounded
    ---------------
    w, h, mdot, T_out (if increasingT), l

    Lower Unbounded
    ---------------
    dQ, T_out (if not increasingT), l

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
            Cf == 0.059/(fluid.rho*v_in*l/fluid.mu)
        ]
