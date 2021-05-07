from math import *
import numpy as np
import matplotlib.pyplot as plt
from ISA import script


C_D_0 = 0.02
C_L = 2.5
e = 0.8
A = 9
h_cruise = 2000

rho = 1.225
V_s = 31.38    #stall speed
n_p = 0.8   #Propellor efficiency

c = 5 #Climb rate
cV = 0.083  #Climb gradient
V = 90 #Cruise speed

sigma=1  #
S_to=750    #Take-off distance
S_l=750     #Landing distance
f=1 #take-off vs landing max weight
power_setting=0.9
cruise_fraction = 1

"""
C_D_0 = 0.0335
C_L = 1.8
e = 0.7
rho = 1.225
V_s = 31.38    #stall speed
n_p = 0.8   #Propellor efficiency
A = 9
c = 5 #Climb rate 
cV = 0.083  #Climb gradient
V = 70 #Cruise speed

sigma= 1  #
S_to = 700    #Take-off distance
S_l= 550     #Landing distance
f=0.95 #take-off vs landing max weight
power_setting=0.9
cruise_fraction = 0.8
"""

def dragcoef(C_D_0, C_L, A, e):
    C_D = C_D_0 + (C_L ** 2) / (pi * A * e)
    return C_D


def Stallload(rho, C_L_max, V_s):
    WS = 0.5 * rho * V_s ** 2 * C_L_max
    return (WS)

def Landing(S_l,rho,C_L_max,f):   #C_L_max for landing, f is fraction between takeoff and landing max weigth
    WS=(C_L_max*rho*(S_l/0.5915)/(2*f))
    return(WS)

def TOPcalc(S_to):
    a = 0.0577
    b = 8.6726
    c = -S_to
    TOP = (-b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    return TOP

def takeoff(x,sigma,C_L_max,S_to):
    C_L_to=C_L_max / 1.21
    y = []
    for i in range(len(x)):
        y_temp=TOPcalc(S_to) / x[i] * C_L_to * sigma
        y.append(y_temp)
    return(y)


def cruise_perf(x, sigma, n_p, power_setting, cruise_fraction, C_D_0, rho, V, A, e, h_cruise):
    sigma = script(h_cruise)/script(0)

    y = []
    for i in range(len(x)):
        y_temp = power_setting / cruise_fraction * n_p * sigma ** 0.75 * (
                (C_D_0 * 0.5 * rho * V ** 3) / (cruise_fraction * x[i]) + (cruise_fraction * x[i]) * 1 / (
                    pi * A * e * 0.5 * rho * V)) ** (-1)
        y.append(y_temp)
    return (y)


def climb_rate(x, n_p, c, A, e, C_D_0, rho):
    y = []
    for i in range(len(x)):
        y_temp = n_p / (c + ((sqrt(x[i]) * sqrt(2 / rho)) / (1.345 * (A * e) ** 0.75 / C_D_0 ** 0.25)))
        y.append(y_temp)
    return (y)


def climb_gradient(x, n_p, cV, C_D_0, C_L_m, rho, A, e):
    C_L = C_L_m / 1.1  #Safety margin of 10% on C_L
    C_D = dragcoef(C_D_0, C_L, A, e)
    y = []
    for i in range(len(x)):
        y_temp = n_p / (sqrt(x[i]) * (cV + C_D/C_L) * sqrt(2 / rho / C_L))
        y.append(y_temp)
    return (y)



x = np.arange(1,3000)
#y = []
#for i in range(len(x)):
#    y_temp = climb_gradient(x[i],0.9,10,0.1,5,1.225, 10, 0.8)
#    y.append(y_temp)



plt.vlines(Stallload(rho, C_L, V_s),0,1, label = "Stall load",color='red')
plt.vlines(Landing(S_l,rho,C_L,f),0,1, label = "Landing")
plt.plot(x,takeoff(x,sigma,C_L,S_to), label = "Take-off - C_L = " + str(C_L))
plt.plot(x,takeoff(x,sigma,C_L-0.3,S_to), label = "Take-off - C_L = " + str(C_L-0.3))
plt.plot(x,takeoff(x,sigma,C_L-0.5,S_to), label = "Take-off - C_L = " + str(C_L-0.5))

plt.plot(x,cruise_perf(x, sigma, n_p, power_setting, cruise_fraction, C_D_0, rho, V, A-3, e, h_cruise), label = "Cruise, A = " + str(A-3))
plt.plot(x,climb_rate(x, n_p, c, A-3, e, C_D_0, rho),label = 'Climb rate, A = ' + str(A-3))
plt.plot(x,climb_gradient(x, n_p, cV, C_D_0, C_L, rho, A-3, e), label='Climb gradient, A = ' + str(A-3))

plt.plot(x,cruise_perf(x, sigma, n_p, power_setting, cruise_fraction, C_D_0, rho, V, A, e, h_cruise), label = "Cruise, A = " + str(A))
plt.plot(x,climb_rate(x, n_p, c, A, e, C_D_0, rho),label = 'Climb rate, A = ' + str(A))
plt.plot(x,climb_gradient(x, n_p, cV, C_D_0, C_L, rho, A, e), label='Climb gradient, A = ' + str(A))

plt.plot(x,cruise_perf(x, sigma, n_p, power_setting, cruise_fraction, C_D_0, rho, V, A+3, e, h_cruise), label = "Cruise, A = " + str(A+3))
plt.plot(x,climb_rate(x, n_p, c, A+3, e, C_D_0, rho),label = 'Climb rate, A = ' + str(A+3))
plt.plot(x,climb_gradient(x, n_p, cV, C_D_0, C_L, rho, A+3, e), label='Climb gradient, A = ' + str(A+3))


plt.ylim(0,0.4)
plt.xlim(0,3000)
plt.legend()
plt.show()
