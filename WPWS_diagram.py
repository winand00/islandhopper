from math import *
import numpy as np
import matplotlib.pyplot as plt
from ISA import script

######
#Write the name of your design option
#Put it then below at the other functions
#Change values and plot
######

class flyingwing:
    def __init__(self):
        # Change these 5 as you wish
        self.C_D_0 = 0.039
        self.C_L = 2.5
        self.e = 0.8
        self.A = 9
        self.h_cruise = 2000

        self.rho = 1.225
        self.V_s = 31.38  # stall speed
        self.n_p = 0.8  # Propellor efficiency

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1

class hydrogen:
    def __init__(self):
        # Change these 5 as you wish
        self.C_D_0 = 0.039
        self.C_L = 2.5
        self.e = 0.75
        self.A = 9
        self.h_cruise = 2000

        self.rho = 1.225
        self.V_s = 31.38  # stall speed
        self.n_p = 0.8  # Propellor efficiency

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1

class claimthisname2:
    def __init__(self):
        # Change these 5 as you wish
        self.C_D_0 = 0.02
        self.C_L = 2.5
        self.e = 0.8
        self.A = 9
        self.h_cruise = 2000

        self.rho = 1.225
        self.V_s = 31.38  # stall speed
        self.n_p = 0.8  # Propellor efficiency

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1

class claimthisname3:
    def __init__(self):
        # Change these 5 as you wish
        self.C_D_0 = 0.02
        self.C_L = 2.5
        self.e = 0.8
        self.A = 9
        self.h_cruise = 2000

        self.rho = 1.225
        self.V_s = 31.38  # stall speed
        self.n_p = 0.8  # Propellor efficiency

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1

class claimthisname4:
    def __init__(self):
        # Change these 5 as you wish
        self.C_D_0 = 0.02
        self.C_L = 2.5
        self.e = 0.8
        self.A = 9
        self.h_cruise = 2000

        self.rho = 1.225
        self.V_s = 31.38  # stall speed
        self.n_p = 0.8  # Propellor efficiency

        self.c = 5  # Climb rate
        self.cV = 0.083  # Climb gradient
        self.V = 90  # Cruise speed

        self.sigma = 1  #
        self.S_to = 750  # Take-off distance
        self.S_l = 750  # Landing distance
        self.f = 1  # take-off vs landing max weight
        self.power_setting = 0.9
        self.cruise_fraction = 1


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
    WS = 0.5 * a.rho * a.V_s ** 2 * CL_value
    return (WS)


def Landing(a, A_value=-1, CL_value=-1):  # C_L_max for landing, f is fraction between takeoff and landing max weigth
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L
    WS = (CL_value * a.rho * (a.S_l / 0.5915) / (2 * a.f))
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
        y_temp = a.n_p / (a.c + ((sqrt(x[i]) * sqrt(2 / a.rho)) / (1.345 * (A_value * a.e) ** 0.75 / a.C_D_0 ** 0.25)))
        y.append(y_temp)
    return (y)


def climb_gradient(a, A_value=-1, CL_value=-1):
    if A_value == -1:
        A_value = a.A
    if CL_value == -1:
        CL_value = a.C_L

    C_L = CL_value / 1.1  # Safety margin of 10% on C_L
    C_D = dragcoef(a, A_value, CL_value)
    y = []
    for i in range(len(x)):
        y_temp = a.n_p / (sqrt(x[i]) * (a.cV + C_D / CL_value) * sqrt(2 / a.rho / CL_value))
        y.append(y_temp)
    return (y)


x = np.arange(1, 3000)
#color = ['firebrick','red','tomato','orange','goldenrod','gold','limegreen','lime','seagreen','violet','magenta','deeppink']

def wpws_plot(a,option = -1):

    #Turn of the options you dont need.

    # Stall load and landing constraints
    plt.vlines(Stallload(a), 0, 1, label="Stall load", color='dimgrey')
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


wpws_plot(flyingwing())

#wpws_plot(claimthisname1())

#wpws_plot(claimthisname2())

#wpws_plot(claimthisname3())

#wpws_plot(claimthisname4())

