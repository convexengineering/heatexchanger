from gpkit import Model, parse_variables, SignomialEquality, SignomialsEnabled, units


class RectangularPipe(Model):
    """
    Defines heat exchanger pipe elements, and the fluid-wall interactions

    Variables
    ---------
    mdot                  [kg/s]   mass flow rate
    w                     [cm]      pipe width
    T_in                  [K]      input temperature
    v_in                  [m/s]    input velocity
    v_out                 [m/s]    output velocity
    P_in         10100    [Pa]     input static pressure
    P_out        10100    [Pa]     output static pressure
    D                     [N]      total drag
    eta_h                 [-]      effectiveness
    eta_h_ref   0.917     [-]      reference effectiveness
    Pf                    [-]      pressure drop parameter
    Pf_ref      21.66     [-]      reference pressure drop parameter
    Re_ref      90550     [-]      reference Reynolds number
    Pr                    [-]      Prandtl number
    fr                    [Pa]     force per frontal area

    Variables of length Nsegments+1
    -------------------------------
    v                     [m/s]    fluid velocity
    T                     [K]      fluid temperature
    P0                    [Pa]     fluid total pressure
    A                     [cm^2]   area

    Variables of length Nsegments
    -----------------------------
    dT                    [K]       Change in fluid temperature over segment
    dQ                    [W]       Magnitude of heat transfer over segment
    T_avg                 [K]       Average temperature over segment
    v_avg                 [m/s]     Average fluid velocity over segment
    l                     [cm]       Reference flow length
    dh                    [cm]       hydraulic diameter
    V_seg                 [cm^3]     Segment volume
    A_seg                 [cm^2]     Segment frontal area
    h_seg                 [cm]       Segment height
    l_seg                 [cm]       Segment flow length
    w_fluid               [cm]      fluid width
    Nu                    [-]       Nusselt number
    Re                    [-]       Reynolds number
    dP                    [Pa]      segment pressure drop
    Tr_int                [K]       wall-fluid interface temperature
    h                     [W/K/m^2] convective heat transfer coefficient
    Cf                    [-]       coefficient of friction

    Upper Unbounded
    --------------
    w, dh, l_seg, V_seg, D, h_seg
    Tr_int (if increasingT), T_in (if not increasingT)

    Lower Unbounded
    ---------------
    w, Re_notlast, dQ, dP
    Tr_int (if not increasingT), T_in (if increasingT)

    """

    def setup(self, Nsegments, Nfins, fluid, increasingT):
        self.fluid = fluid
        self.increasingT = increasingT

        exec parse_variables(RectangularPipe.__doc__)
        self.Re_notlast = Re[:-1]  # unbounded b/c only Re[-1] is fit

        temp = [T_avg**2 == T[1:] * T[:-1],
                T_in == T[0]]

        # TODO: why is Tr_int not just equal to T[1:]?

        if increasingT:
            temp.extend([T[1:] >= T[:-1] + dT,
                         # definition of effectiveness
                         dT * eta_h**-1 + T[:-1] <= Tr_int,
                         Tr_int >= T[1:]])
        else:
            temp.extend([T[:-1] >= T[1:] + dT,
                         dT * eta_h**-1 + Tr_int <= T[:-1],
                         Tr_int <= T[1:]])

        with SignomialsEnabled():
            # except where noted these constraints become posynomial
            flow = [
                mdot == Nfins * fluid.rho * v * A,  # mass conservation
                v_in == v[0],
                v_out == v[-1],
                v_avg**2 == v[:-1] * v[1:],
                # force per frontal area
                fr == Pf * (0.5 * fluid.rho * v_in**2),
                # inlet total pressure # SIGNOMIAL
                P0[0] <= P_in + 0.5 * fluid.rho * v_in**2,
                # exit total pressure
                P0[-1] >= P_out + 0.5 * fluid.rho * v_out**2,
                P0[0] >= P0[-1] + 0.5 * fluid.rho * v_in**2 * Pf,
                P0[:-1] >= P0[1:] + dP,
                dP <= mdot * (v[:-1]/A[:-1] - v[1:]/A[1:]),
                A_seg**2 == A[:-1]*A[1:],
                dP == 0.5 * fluid.rho * v_avg**2 * Cf * l_seg / dh,

                # effectiveness fit
                eta_h / eta_h_ref == 0.799 * (Re[-1]/Re_ref)**-0.0296,
                eta_h <= 0.844,  # boundary to make sure fit is valid

                # pressure drop fit
                (Pf/Pf_ref)**0.155 >= (0.475*(Re[-1]/Re_ref)**0.00121
                                       + 0.0338*(Re[-1]/Re_ref)**-0.336),

                D >= fr * Nfins * A[0]
            ]

        # Geometry definitions
        geom = [A_seg == w_fluid * h_seg,  # cross sectional area of channel
                V_seg == Nfins * A_seg * l_seg,  # total volume of all channels
                # hydraulic diameter with geometric mean approximation
                dh * (w_fluid * h_seg)**0.5 == 2 * A_seg,
                h_seg >= 0.1*units('cm'),
                ]
        with SignomialsEnabled():
            for i in range(Nsegments):
                geom.extend([l[i] <= l_seg[:i + 1].sum()])

        # Friction and heat transfer
        friction = [dQ <= mdot * fluid.c * dT,
                    Cf**5 * Re == (0.059)**5,
                    Re == fluid.rho * v_avg * l / fluid.mu,
                    Pr == fluid.mu * fluid.c / fluid.k,
                    # Nusselt number definition (fully turbulent)
                    Nu == 0.0296 * Re**(4./5) * Pr**(1./3),
                    h * l == Nu * fluid.k,
                    ]

        return [fluid, temp, flow, geom, friction]
