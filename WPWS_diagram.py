from math import *
import numpy as np
import matplotlib.pyplot as plt
from tool import tool


# from ISA import script


def script(h):
    # ISA CALCULATOR WINAND E4G#
    #############INPUT

    '''
    print("     **** ISA calculator ****")
    print()
    print("1. Calculate ISA for altitude in meters")
    print("2. Calculate ISA for altitude in feet")
    print("3. Calculate ISA for altitude in FL")
    print()
    c = float(input("Enter your choice:"))

    #############CHOICES

    if c == 1.:
        print()
        h = float(input("Enter altitude [m]: "))

    elif c == 2.:
        print()
        h = float(input("Enter altitude [ft]: "))
        h = h * 0.3048

    elif c == 3.:
        print()
        h = float(input("Enter altitude [FL]: "))
        h = h * 30.48

    #############TEMPERATURE

    print()
    print("                        Sea level temperature")
    t0 = (input("Press enter when using standard temperature):"))
    if t0 == "":
        t0 = 288.15
    else:
        t0 = float(t0)
    '''

    t0 = 288.15
    #############CONSTANTS

    h1 = min(h, 11000.)
    h2 = min(h, 20000.)
    h3 = min(h, 32000.)
    h4 = min(h, 47000.)
    h5 = min(h, 51000.)
    h6 = min(h, 71000.)
    h7 = min(h, 84852.)

    g0 = 9.80665
    r = 287.05
    p0 = 101325
    rho0 = 1.225
    a1 = -0.0065
    a2 = 0.0000
    a3 = 0.0010
    a4 = 0.0028
    a5 = 0.0000
    a6 = -0.0028
    a7 = -0.0020

    #############VALUES
    t1 = t0 + h1 * a1
    p1 = p0 * (t1 / t0) ** (-g0 / (a1 * r))
    rho1 = rho0 * (t1 / t0) ** (-g0 / (a1 * r) - 1)

    t2 = t1 + (h2 - h1) * a2
    p2 = p1 * e ** (-g0 / (r * t2) * (h2 - h1))
    rho2 = rho1 * e ** (-g0 / (r * t2) * (h2 - h1))

    t3 = t2 + (h3 - h2) * a3
    p3 = p2 * (t3 / t2) ** (-g0 / (a3 * r))
    rho3 = rho2 * (t3 / t2) ** (-g0 / (a3 * r) - 1)

    t4 = t3 + (h4 - h3) * a4
    p4 = p3 * (t4 / t3) ** (-g0 / (a4 * r))
    rho4 = rho2 * (t4 / t3) ** (-g0 / (a4 * r) - 1)

    t5 = t4 + (h5 - h4) * a5
    p5 = p4 * e ** (-g0 / (r * t5) * (h5 - h4))
    rho5 = rho4 * e ** (-g0 / (r * t5) * (h5 - h4))

    t6 = t5 + (h6 - h5) * a6
    p6 = p5 * (t6 / t5) ** (-g0 / (a6 * r))
    rho6 = rho5 * (t6 / t5) ** (-g0 / (a6 * r) - 1)

    t7 = t6 + (h7 - h6) * a7
    p7 = p6 * (t7 / t6) ** (-g0 / (a7 * r))
    rho7 = rho6 * (t7 / t6) ** (-g0 / (a7 * r) - 1)

    #############

    tt = [t0, t1, t2, t3, t4, t5, t6, t7]
    pp = [p0, p1, p2, p3, p4, p5, p6, p7]
    rr = [rho0, rho1, rho2, rho3, rho4, rho5, rho6, rho7]

    #############sorry

    if h > 84852:
        print()
        print("not yet possible :( ")

    #############CALCULATING

    else:
        if h <= 11000:
            oo = 1
        elif 11000 < h <= 20000:
            oo = 2
        elif 20000 < h <= 32000:
            oo = 3
        elif 32000 < h <= 47000:
            oo = 4
        elif 47000 < h <= 51000:
            oo = 5
        elif 51000 < h <= 71000:
            oo = 6
        elif 71000 < h <= 84852:
            oo = 7
        # print()
        tc = tt[oo] - 273.15
        # print("Temperature : ", round(tt[oo], 2), "K", " (", round(tc, 2), "'C)")
        ps = (pp[oo] / p0) * 100
        # print("Pressure    : ", round(pp[oo], 2), "Pa", "(", round(ps, 2), "%SL)")
        rhos = (rr[oo] / rho0) * 100
        # print("Density     : ", round(rr[oo], 2), "kg/m3", "(", round(rhos, 2), "%SL)")\
    return (rr[oo])


######
# Write the name of your design option
# Put it then below at the other functions
# Change values and plot
######

class flyingwing:
    def __init__(self):
        # Change these first 7 as you wish
        self.C_D_0 = 0.007
        self.C_L = 2.0
        self.C_L_cruise = 0.8
        self.e = 0.8
        self.A = 7

        self.h_cruise = 10000*0.3048
        self.m_energy = 2495.1275 #[kg]

        self.battery = True  # Aircraft on batteries

        # Change these engine variables as you wish
        self.D = 0 #Propellor diameter
        self.B = 4 #number of blader per propellor
        self.N = 4 #Number of engines
        self.efficiency_r = 0.1 #Fraction of total used energy that is recovered for other systems

        self.rho0 = 1.225
        self.rho= script(self.h_cruise)
        self.V_s = 43  # stall speed
        self.n_p = 0.8  # Propellor efficiency
        self.C_L_takeoff = self.C_L/(1.1**2)

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = self.S_to  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1
        self.W = 8618.255 * 9.80665  # kg

        self.S = 0  #Wing surface area


        self.specific_energy = 550*3600 #Specific energy of fuel [J/kg]
        self.efficiency_fuelcell = 0.9   # Efficiency fuel cell

        self.P = 0 #Max power [W]
        self.L_over_D= self.C_L_cruise/dragcoef(self,CL_value=self.C_L_cruise)




class hydrogen:
    def __init__(self):
        # Change these 7 as you wish
        self.C_D_0 = 0.039
        self.C_L = 2.2
        self.C_L_cruise = 0.8
        self.e = 0.8
        self.A = 9
        self.h_cruise = 3048
        self.m_energy = 104.8 #[kg]
        self.battery = False  # Aircraft on batteries

        # Change these engine variables as you wish
        self.D = 0  # Propellor diameter
        self.B = 4  # number of blader per propellor
        self.N = 2  # Number of engines
        self.efficiency_r = 0.1  # Fraction of total used energy that is recovered for other systems

        self.rho0 = 1.225
        self.rho = script(self.h_cruise)
        self.V_s = 43  # stall speed
        self.n_p = 0.85  # Propellor efficiency
        self.C_L_takeoff = self.C_L / (1.1 ** 2)

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1

        self.W = 8618.255*9.80655  # N


        self.S = 0  # Wing surface area
        self.specific_energy = 120000000  # Specific energy of fuel [J/kg]
        self.efficiency_fuelcell = 0.6   # Efficiency fuel cell
        self.P = 0  # Max power [W]
        self.L_over_D = self.C_L_cruise / dragcoef(self,CL_value=self.C_L_cruise)

class conc_batteries:
    def __init__(self):
        self.C_D_0 = 0.025
        self.C_L = 2.4
        self.C_L_cruise = 0.8
        self.e = 0.85
        self.A = 12
        self.h_cruise = 3048
        self.m_energy = 1360  # [kg]
        self.battery = True  # Aircraft on batteries

        # Change these engine variables as you wish
        self.D = 2.77  # Propeller diameter
        self.B = 6 # number of blades per propeller
        self.N = 2  # Number of engines
        self.efficiency_r = 0.1  # Fraction of total used energy that is recovered for other systems

        self.rho0 = 1.225
        self.rho = script(self.h_cruise)
        self.V_s = 43  # stall speed
        self.n_p = 0.8  # Propellor efficiency
        self.C_L_takeoff = self.C_L / (1.1 ** 2)

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.5
        self.cruise_fraction = 1
        self.W = 8618.255*9.80655  # N

        self.S = 0  # Wing surface area
        self.specific_energy = 600*3600  # Specific energy of fuel [J/kg]
        self.efficiency_fuelcell = 1   # Efficiency fuel cell
        self.P = 0  # Max power [W]
        self.L_over_D = self.C_L_cruise / dragcoef(self,CL_value=self.C_L_cruise)


class distributed:
    def __init__(self):
        self.C_D_0 = 0.0307
        self.C_L = 3.66
        self.C_L_cruise = 1.08
        self.e = 0.8
        self.A = 15.16
        self.h_cruise = 3000
        self.m_energy = 2000  # [kg]
        self.battery = True  # Aircraft on batteries

        # Change these engine variables as you wish
        self.D = 1.72  # Propellor diameter
        self.B = 4  # number of blader per propellor
        self.N = 12  # Number of engines
        self.efficiency_r = 0.1  # Fraction of total used energy that is recovered for other systems

        self.rho0 = 1.225
        self.rho = script(self.h_cruise)
        self.V_s = 84  # stall speed
        self.n_p = 0.8  # Propellor efficiency
        self.C_L_takeoff = self.C_L / (1.1 ** 2)

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1
        self.W = 8618.255*9.80655  # N

        self.S = 29.74  # Wing surface area
        self.specific_energy = 550*3600  # Specific energy of fuel [J/kg]
        self.efficiency_fuelcell = 1   # Efficiency fuel cell

        self.P = 1293925  # Max power [W]
        self.L_over_D = self.C_L_cruise / dragcoef(self,CL_value=self.C_L_cruise)


class claimthisname4:
    def __init__(self):
        # Change these 5 as you wish
        self.C_D_0 = 0.02
        self.C_L = 2.0
        self.e = 0.8
        self.A = 7
        self.h_cruise = 10000*0.3048
        self.m_energy = 2000  # [kg]
        self.battery = False  # Aircraft on batteries

        # Change these engine variables as you wish
        self.D = 2.69  # Propellor diameter
        self.B = 4  # number of blader per propellor
        self.N = 2  # Number of engines
        self.efficiency_r = 0.1  # Fraction of total used energy that is recovered for other systems

        self.rho0 = 1.225
        self.rho = script(self.h_cruise)
        self.V_s = 31.38  # stall speed
        self.n_p = 0.8  # Propellor efficiency
        self.C_L_takeoff = self.C_L / (1.1 ** 2)

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1
        self.W = 8618.255*9.80655  # N

        self.S = 0  # Wing surface area
        self.specific_energy = 46200000  # Specific energy of fuel [J/kg]
        self.efficiency_fuelcell = 0.9   # Efficiency fuel cell
        self.P = 0  # Max power [W]
        self.L_over_D = self.C_L / dragcoef(self,CL_value=self.C_L_cruise)


def dragcoef(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L
    C_D = a.C_D_0 + (CL_value ** 2) / (pi * A_value * a.e)
    return C_D


def Stallload(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L
    WS = 0.5 * a.rho0 * a.V_s ** 2 * CL_value
    return (WS)


def Landing(a, A_value=-1, CL_value=-1):  # C_L_max for landing, f is fraction between takeoff and landing max weigth
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L
    WS = (CL_value * a.rho0 * (a.S_l / 0.5915) / (2 * a.f))
    return (WS)


def TOPcalc(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L
    A = 0.0577
    B = 8.6726
    C = -(a.S_to)
    TOP = (-B + sqrt(B ** 2 - 4 * A * C)) / (2 * A)
    return TOP


def takeoff(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L
    C_L_to = CL_value / 1.21
    y = []
    for i in range(len(x)):
        y_temp = TOPcalc(a, A_value, CL_value) / x[i] * C_L_to * a.sigma
        y.append(y_temp)
    return (y)


def cruise_perf(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L

    sigma = script(a.h_cruise) / script(0)

    y = []
    for i in range(len(x)):
        y_temp = a.power_setting / a.cruise_fraction * a.n_p * sigma ** 0.75 * (
                (a.C_D_0 * 0.5 * a.rho * a.V ** 3) / (a.cruise_fraction * x[i]) + (a.cruise_fraction * x[i]) * 1 / (
                pi * A_value * a.e * 0.5 * a.rho * a.V)) ** (-1)
        y.append(y_temp)
    return (y)


def climb_rate(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L

    y = []
    for i in range(len(x)):
        y_temp = a.n_p / (a.c + ((sqrt(x[i]) * sqrt(2 / a.rho0)) / (1.345 * (A_value * a.e) ** 0.75 / a.C_D_0 ** 0.25)))
        y.append(y_temp)
    return (y)


def climb_gradient(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L


    C_L = sqrt(3 * a.C_D_0 * pi * a.A * a.e)
    C_D = dragcoef(a, A_value, CL_value)
    y = []
    for i in range(len(x)):
        y_temp = a.n_p / (sqrt(x[i]) * (a.cV + C_D / CL_value) * sqrt(2 / a.rho0 / CL_value))
        y.append(y_temp)
    return (y)


x = np.arange(1, 3000)


# color = ['firebrick','red','tomato','orange','goldenrod','gold','limegreen','lime','seagreen','violet','magenta','deeppink']

def wpws_plot(a, option=-1):
    # Turn of the options you dont need.

    # Stall load and landing constraints
    plt.vlines(Stallload(a), 0, 1, label="Stall load", color='red')
    plt.vlines(Landing(a), 0, 1, label="Landing", color='dimgray')

    # Take-off constraints, varying CL_max
    plt.plot(x, takeoff(a), label="Take-off - C_L = " + str(a.C_L), color='firebrick')
    plt.plot(x, takeoff(a, CL_value=a.C_L - 0.3), label="Take-off - C_L = " + str(a.C_L - 0.3), color='red')
    plt.plot(x, takeoff(a, CL_value=a.C_L - 0.5), label="Take-off - C_L = " + str(a.C_L - 0.5), color='tomato')

    # Cruise constraints, varying A
    plt.plot(x, cruise_perf(a, A_value=a.A - 3), label="Cruise, A = " + str(a.A - 3), color='limegreen')
    plt.plot(x, cruise_perf(a), label="Cruise, A = " + str(a.A), color='cornflowerblue')
    plt.plot(x, cruise_perf(a, A_value=a.A + 3), label="Cruise, A = " + str(a.A + 3), color='blueviolet')

    # Climb rate constraints, Varying A
    plt.plot(x, climb_rate(a, A_value=a.A - 3), label='Climb rate, A = ' + str(a.A - 3), color='lime')
    plt.plot(x, climb_rate(a), label='Climb rate, A = ' + str(a.A), color='cyan')
    plt.plot(x, climb_rate(a, A_value=a.A + 3), label='Climb rate, A = ' + str(a.A + 3), color='magenta')

    # Climb gradient constraints, Varying A
    plt.plot(x, climb_gradient(a, A_value=a.A - 3), label='Climb gradient, A = ' + str(a.A - 3), color='yellowgreen')
    plt.plot(x, climb_gradient(a), label='Climb gradient, A = ' + str(a.A), color='dodgerblue')
    plt.plot(x, climb_gradient(a, A_value=a.A + 3), label='Climb gradient, A = ' + str(a.A + 3), color='deeppink')

    plt.ylim(0, 0.4)
    plt.xlim(0, 3000)
    plt.legend()
    plt.show()


def design_point(a, WS, WP):
    print('Density at cruise altitude is =',script(a.h_cruise))
    S_l = WS * (2 * a.f) / (a.C_L * a.rho0) * 0.5915
    Power = a.W / WP

    a.P = Power
    a.S = a.W / WS
    print('Surface area =', a.S)
    print('Landing length = ', S_l)
    print('Power required = ', Power)
    a.D = 0.55*((Power/(1000*a.N))**(0.25))
    print('Propeller diamter in [m] = ',a.D,'With ',a.N,' engines')
    sigma = script(a.h_cruise) / script(0)

    power_setting = ((a.C_D_0 * 0.5 * a.rho * a.V ** 3) / (a.cruise_fraction * WS) + (a.cruise_fraction * WS) * 1 / (
                pi * a.A * a.e * 0.5 * a.rho * a.V)) * WP * a.cruise_fraction * (a.n_p * sigma ** 0.75) ** (-1)
    cruisepower = Power * power_setting
    print('Cruise power = ', cruisepower)

    c = a.n_p / WP - (((sqrt(WS) * sqrt(2 / a.rho0)) / (1.345 * (a.A * a.e) ** 0.75 / a.C_D_0 ** 0.25)))
    print('Climb rate = ', c)


    C_L = sqrt(3 * a.C_D_0 * pi * a.A * a.e)
    C_D = dragcoef(a, CL_value=C_L)
    cV = a.n_p * (1 / WP) * (1 / (sqrt(WS * 2 / a.rho0 / C_L))) - C_D / C_L
    print('Climb gradient = ', degrees(atan(cV)))
    print('Lift over drag is ',a.L_over_D)

def Tool(a,WS, WP):
    plt.scatter(WS,WP)
    wpws_plot(a)
    design_point(a,WS, WP)
    tool(a.C_D_0, a.A, a.e, a.W, a.rho, a.rho0, a.S, a.specific_energy, a.m_energy, a.W / 9.80665, a.L_over_D,
         a.efficiency_fuelcell, a.n_p, a.P, a.C_L_takeoff, a.C_L, a.D, a.B, a.N, a.efficiency_r, a.battery)




####
#Fill in aircraftname, WS, WP
#WS = 1553
#WP = 0.0653
#Tool(hydrogen(),WS,WP)
#wpws_plot(hydrogen())

WS = 1408.6
WP = 0.07215
Tool(flyingwing(),WS,WP)




WS = 2842.4
WP = 0.06535
Tool(distributed(),WS,WP)




