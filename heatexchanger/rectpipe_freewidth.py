from gpkit import Model, parse_variables


class RectangularPipe(Model):
    """
    Variables
    ---------
    mdot     [kg/s]   mass flow rate
    T_in     [K]      input temperature
    v_in     [m/s]    input velocity
    h        [m]      pipe height

    Variables of length Nsegments+1
    -------------------------------
    v        [m/s]    fluid velocity
    T        [K]      fluid temperature
    w        [m]      pipe width
    A        [m^2]    pipe area

    Variables of length Nsegments
    -----------------------------
    dT       [K]      Change in fluid temperature over segment
    dQ       [W]      Magnitude of heat transfer over segment

    Upper Unbounded
    ---------------
    v_after_input, w, h, mdot, T_out (if increasingT)

    Lower Unbounded
    ---------------
    w, dQ, T_out (if not increasingT)

    """
    def setup(self, Nsegments, fluid, increasingT):
        exec parse_variables(RectangularPipe.__doc__)
        self.increasingT = increasingT
        self.T_out = T[-1]
        self.v_after_input = v[1:]
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
        ]
