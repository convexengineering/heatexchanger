from gpkit import Model, parse_variables
import numpy as np


class Channel(Model):
    """Single channel HX model

    Variables
    ---------
    A_r                  [m^2]    frontal area
    h                    [m]      channel height
    l                    [m]      channel length
    eta_h_ref      0.917 [-]      reference effectiveness
    Pf_ref        21.66  [-]      reference pressure drop parameter
    Re_ref     90550     [-]      reference Reynolds number

    Upper Unbounded
    ---------------
    A_r, l

    Lower Unbounded
    ---------------
    A_r, l

    """
    def setup(self):
        exec parse_variables(Channel.__doc__)
        return [A_r == h**2]

    def dynamic(self, state):
        return ChannelP(self, state)


class ChannelP(Model):
    """Single channel air HX performance model

    Variables
    ---------
    cp       1004   [J/(kg*K)]    heat capacity of air
    dT              [-]           wall/free-stream temp ratio - 1
    eta_h           [-]           effectiveness
    fV              [-]           air velocity ratio across channel
    fr              [N/m^2]       force per frontal area
    Hdot        9.9 [J/s]         heat flow rate
    mdot            [kg/s]        air mass flow rate
    Pf              [-]           pressure drop parameter
    Re              [-]           Reynolds number
    Tr     85+273   [K]           wall temperature
    A_e             [m^2]         flow exit area


    Upper Unbounded
    ---------------
    fV, A_e, P0_i, mu_o, mu_i

    Lower Unbounded
    ---------------
    fr, fV, V_o, P_o, V_i, rho_i, rho_o, mu_o, mu_i, Pf_ref

    """
    def setup(self, channel, state):
        exec parse_variables(ChannelP.__doc__)
        self.channel = channel

        mu_i = self.mu_i = state.mu_i
        mu_o = self.mu_o = state.mu_o
        rho_i = self.rho_i = state.rho_i
        rho_o = self.rho_o = state.rho_o
        V_i = self.V_i = state.V_i
        V_o = self.V_o = state.V_o
        T_i = self.T_i = state.T_i
        T_o = self.T_o = state.T_o
        P0_i = self.P0_i = state.P0_i
        P_o = self.P_o = state.P_o
        A_r = self.channel.A_r
        eta_h_ref = self.eta_h_ref = self.channel.eta_h_ref
        Pf_ref = self.Pf_ref = self.channel.Pf_ref
        Re_ref = self.channel.Re_ref
        l = self.l = self.channel.l

        Pf_rat = Pf/Pf_ref
        Re_rat = Re/Re_ref
        eta_h_rat = eta_h/eta_h_ref
        alpha = T_o/T_i

        return [mdot == rho_i*V_i*A_r,
                mdot == rho_o*V_o*A_e,
                Hdot == mdot*cp*dT*eta_h*T_i,
                fV == V_o/V_i,
                # channel Reynolds number (geometric average)
                Re == (rho_i*rho_o*V_i*V_o/mu_i/mu_o)**(0.5)*self.channel.l,

                # effectiveness fit
                eta_h/eta_h_ref == 0.799*Re_rat**-0.0296,
                eta_h <= 0.844,  # boundary to make sure fit is valid
                # pressure drop fit
                Pf_rat**0.155 >= 0.475*Re_rat**0.00121 + 0.0338*Re_rat**-0.336,

                alpha >= 1 + dT*eta_h,
                fr == Pf*(0.5*rho_i*V_i**2),  # def'n of fr
                dT*T_i + T_i <= Tr,
                P0_i >= P_o + 0.5*rho_i*V_i**2*Pf + 0.5*rho_o*V_o**2,
                T_o <= Tr,
                ]


class HXState(Model):
    """HX flow state model

    Variables
    ---------
    mu_i      1.7e-5                [kg/(m*s)]  incoming dynamic viscosity
    mu_o      1.7e-5                [kg/(m*s)]  exit dynamic viscosity
    P0_i      0.5*1.25*30**2+101000 [Pa]        incoming total pressure
    P_i                             [Pa]        incoming static pressure
    P_o       90000                 [Pa]        exit static pressure
    R         287.1                 [J/(kg*K)]  specific gas constant of air
    rho_i                           [kg/(m^3)]  incoming air density
    rho_o                           [kg/(m^3)]  exiting air density
    V_i       30                    [m/s]       incoming air velocity
    V_o                             [m/s]       exiting air velocity
    T_i       -20+273               [K]         incoming air temperature
    T_o                             [K]         exiting air temperature

    Upper Unbounded
    ---------------
    V_o, T_o

    Lower Unbounded
    ---------------
    rho_o, rho_i, P_i

    """
    # calc_p0in = lambda self, c: c[self.P_i] + 0.5*c[self.rho_i]*c[self.V_i]**2

    def setup(self):
        exec parse_variables(HXState.__doc__)
        return [P0_i >= P_i + 0.5*rho_i*V_i**2,
                V_o == V_i*T_o/T_i,
                P_o == rho_o*R*T_o,
                P_i == rho_i*R*T_i,
                T_o >= T_i]


class HX(Model):
    """Heat eXchanger model

    Upper Unbounded
    ---------------
    mdot, A_e, A_r

    Lower Unbounded
    ---------------
    fr, rho_i, P_i
    """
    def setup(self, state):
        self.channel = Channel()
        self.channelP = self.channel.dynamic(state)
        self.state = state

        self.mdot = self.channelP.mdot
        self.fr = self.channelP.fr
        self.A_e = self.channelP.A_e
        self.A_r = self.channel.A_r
        self.rho_i = self.state.rho_i
        self.P_i = self.state.P_i

        return self.channel, self.channelP, self.state


if __name__ == "__main__":
    state = HXState()
    m = HX(state)
    m.cost = m.channel.A_r*m.channelP.A_e*m.channel.l*m.channelP.fr
    sol = m.solve()
