from math import *
import numpy as np
import matplotlib.pyplot as plt

# --------- Green Indicators Tool --------- #
# Before using this tool make sure that the parameters used are the latest version of the parameters #

# -- Parameter List -- #

# GI1: Design efficiency #

class Parameters:
    def __init__(self):
        #General Parameters
        self.MTOW = 8618        # Maximum Take-off Weight [kg]

        #Efficiency Parameters
        self.CL = 1.3           #CL [-] (Generally use cruise conditions)
        self.C_D_0 = 0.03       #CD_0 [-]
        self.A = 10             #Aspect ratio [-]
        self.e = 0.8            #Oswald effiency [-]
        self.dW = 1             #Extra added weight [kg] (use negative value for weight reduction)
        self.G = 4.73           #Weight Growth Factor [-] (Get this value from the Iteration_Method_tool.xcl file)
        self.n_prop = 0.95      #Propeller efficiency [-]
        self.n_engine = 0.93    #Engine efficiency [-]
        self.n_pmad = 0.9       #PMDAD efficiency [-]
        self.n_cooling = 0.4    #Cooling efficiency [-]
        self.n_fuelcell = 0.65  #Fuel cell efficiency [-]

        #Noise Parameters
        self.Pto = 1100000      #Take-off Power [W]
        self.Dp = 3.0           #Propeller diameter [m]
        self.B = 5              #Number of blades per propeller [-]
        self.nrpm = 2500        #Engine RPM [evolutions/minute]
        self.c = 343            #Speed of sound [m/s]
        self.N = 2              #Number of propellers [-]
        self.r = 350            #Distance to observer [m]





def Efficiency(a):
    C_D = a.C_D_0 + (a.CL ** 2) / (pi * a.A* a.e)
    GI1 = a.CL / C_D * 1/(a.MTOW + a.dW * a.G) * a.n_prop * a.n_engine * a.n_pmad * a.n_cooling * a.n_fuelcell
    print("Efficiency indicator =", GI1)
    return GI1

def Noise(a):
    P_br = a.P_TO / a.n_prop
    M_t = a.n_p * pi * D_p / (60*a.c)
    SPL_max = 83.4 + 15.3 * log10(P_br/1000) - 20 * log10(D_p) + 38.5 * M_t - 3 * (a.B-2) + 10 * log10(N) - 20 * log(a.r)
    return SPL_max
def Recyclability(a):
    RC = 1
    return RC

design = Parameters1()

g_1 = Efficiency(design)