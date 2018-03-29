from gpkit import Model, parse_variables, SignomialEquality, SignomialsEnabled, units


class RectangularPipe(Model):
    """
    Defines heat exchanger pipe elements, and the fluid-wall interactions

    Variables
    ---------
    mdot                  [kg/s]   mass flow rate
    w                     [m]      pipe width
    T_in                  [K]      input temperature
    v_in                  [m/s]    input velocity
    v_out                 [m/s]    output velocity
    P_in        101000    [Pa]     input static pressure
    P_out       101000    [Pa]     output static pressure
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

    Variables of length Nsegments
    -----------------------------
    dT                    [K]       Change in fluid temperature over segment
    dQ                    [W]       Magnitude of heat transfer over segment
    T_avg                 [K]       Average temperature over segment
    v_avg                 [m/s]     Average fluid velocity over segment
    l                     [m]       Reference flow length
    dh                    [m]       hydraulic diameter
    V_seg                 [m^3]     Segment volume
    A_seg                 [m^2]     Segment frontal area
    h_seg                 [m]       Segment height
    l_seg                 [m]       Segment flow length
    Nu                    [-]       Nusselt number
    Re                    [-]       Reynolds number
    dP                    [Pa]      segment pressure drop
    Tr_int                [K]       wall-fluid interface temperature
    h                     [W/K/m^2] convective heat transfer coefficient
    Cf                    [-]       coefficient of friction

    Upper Unbounded
    --------------
    w, dh, l_seg, V_seg, D, Tr_int (if increasingT), h_seg

    Lower Unbounded
    ---------------
    w, Nu_notlast, dQ, Tr_int (if not increasingT), dP

    """

    def setup(self, Nsegments, fluid, increasingT):
        self.fluid = fluid
        self.increasingT = increasingT

        exec parse_variables(RectangularPipe.__doc__)
        self.Nu_notlast = Nu[:-1]

        temp = [T_avg**2 == T[1:]*T[:-1],
                T_in == T[0]]

        if increasingT:
            temp.extend([T[1:] >= T[:-1] + dT,
                         dT*eta_h**-1 + T[0:-1] <= Tr_int,  # definition of effectiveness
                         Tr_int >= T[1:]])
        else:
            temp.extend([T[:-1] >= T[1:] + dT,
                         dT*eta_h**-1 + Tr_int <= T[0:-1],
                         Tr_int <= T[1:]])

        Pf_rat = Pf/Pf_ref
        Re_rat = Re/Re_ref
        eta_h_rat = eta_h/eta_h_ref
        alpha = T[1:]/T[:-1]

        with SignomialsEnabled():
            flow = [mdot == fluid.rho*v_avg*A_seg, # mass conservation
                    v_in == v[0],
                    v_out == v[-1],
                    v_avg**2 == v[0:-1]*v[1:],
                    fr == Pf*(0.5*fluid.rho*v_in**2),  # force per frontal area
                    P0[0] <= P_in + 0.5*fluid.rho*v_in**2, # inlet total pressure # signomial
                    P0[-1] >= P_out + 0.5*fluid.rho*v_out**2, # exit total pressure
                    P0[0] >= P0[-1] + 0.5*fluid.rho*v_in**2*Pf,
                    P0[:-1] >= P0[1:] + dP,
                    #dP*Nsegments == fr,
                    dP <= fluid.rho*v[0:-1]*(v[0:-1] - v[1:]),
                    dP == 0.5*fluid.rho*v_avg**2*Cf*l_seg/dh,

                    # effectiveness fit
                    eta_h/eta_h_ref == 0.799*Re_rat[-1]**-0.0296,
                    eta_h <= 0.844,  # boundary to make sure fit is valid

                    # pressure drop fit
                    Pf_rat**0.155 >= 0.475*Re_rat[-1]**0.00121 + 0.0338*Re_rat[-1]**-0.336,

                    D >= fr*A_seg[0]
                    ]  # turns into a posynomial

        # Geometry definitions
        geom = [A_seg == w*h_seg,
                V_seg == A_seg*l_seg,
                dh*(w*h_seg)**0.5 == 2*A_seg,   # hydraulic diameter with geometric mean approximation
               ]
        with SignomialsEnabled():
            for i in range(Nsegments):
                geom.extend([l[i] <= l_seg[0:i+1].sum()])

        # Friction and heat transfer
        friction = [dQ  <= mdot*fluid.c*dT,
                    Cf**5*Re == (0.059)**5,
                    Re  == (fluid.rho*v_avg*l/fluid.mu),
                    Pr  == fluid.mu*fluid.c/fluid.k,
                    Nu  == 0.0296*Re**(4./5.)*Pr**(1./3.),   # Nusselt number definition (fully turbulent)
                    h*l == Nu*fluid.k,
                    ]

        return [fluid, temp, flow, geom, friction]
