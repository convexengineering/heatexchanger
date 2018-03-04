from gpkit import Model, parse_variables, SignomialEquality, SignomialsEnabled, units


class RectangularPipe(Model):
    """
    Variables
    ---------
    mdot                  [kg/s]   mass flow rate
    w                     [m]      pipe width
    h                     [m]      pipe height
    dh                    [m]      hydraulic diameter
    A                     [m^2]    pipe area
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
    fr                    [Pa]     force per frontal area
    
    Variables of length Nsegments+1
    -------------------------------
    v                     [m/s]    fluid velocity
    T                     [K]      fluid temperature
    P0                    [Pa]     fluid total pressure

    Variables of length Nsegments
    -----------------------------
    dT                    [K]      Change in fluid temperature over segment
    dQ                    [W]      Magnitude of heat transfer over segment
    v_avg                 [m/s]    Average fluid velocity over segment
    l                     [m]      Reference flow length
    Cf                    [-]      Coefficient of friction over segment 
    Re                    [-]      Reynolds number
    D_seg                 [N]      segment drag
    dP                    [Pa]     segment pressure drop
  
    Upper Unbounded
    ---------------
    w, h, mdot, T_out (if increasingT), l, P0

    Lower Unbounded
    ---------------
    dQ, T_out (if not increasingT), P0

    """
    def setup(self, Nsegments, fluid, increasingT):
        exec parse_variables(RectangularPipe.__doc__)
        self.increasingT = increasingT
        self.T_out = T[-1]
        if increasingT:
            temp = T[1:] >= T[:-1] + dT
        else:
            temp = T[:-1] >= T[1:] + dT

        Pf_rat = Pf/Pf_ref
        Re_rat = Re/Re_ref
        eta_h_rat = eta_h/eta_h_ref
        alpha = T[1:]/T[:-1]

        with SignomialsEnabled():
            pressure = [#alpha >= 1 + dT*eta_h,
                        v_in == v[0],
                        v_out == v[-1],
                        fr == Pf*(0.5*fluid.rho*v_in**2),  # force per frontal area
                        P0[0] <= P_in + 0.5*fluid.rho*v_in**2, # inlet total pressure
                        P0[-1] >= P_out + 0.5*fluid.rho*v_out**2, # exit total pressure 
                        P0[0] >= P0[-1] + 0.5*fluid.rho*v_out**2*Pf + 0.5*fluid.rho*v_in**2,
                        P0[:-1] >= P0[1:] + dP,
                        # effectiveness fit
                        eta_h/eta_h_ref == 0.799*Re_rat[-1]**-0.0296,
                        eta_h <= 0.844,  # boundary to make sure fit is valid
                        # pressure drop fit
                        Pf_rat**0.155 >= 0.475*Re_rat[-1]**0.00121 + 0.0338*Re_rat[-1]**-0.336,


                        ]

        return [
            fluid, temp, 
            pressure,
            T[0] == T_in, v[0] == v_in,
            mdot == fluid.rho*v*A,
            A == w*h,
            dh*(w*h)**0.5 == 2*A, # hydraulic diameter with geometric mean approximation
            dQ <= mdot*fluid.c*dT,
            v_avg**2 == v[0:-1]*v[1:],
            Re == (fluid.rho*v_avg*l/fluid.mu),
            Cf**5*Re == (0.059)**5,
        ]
