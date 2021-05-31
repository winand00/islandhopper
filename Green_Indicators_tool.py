from math import *
import numpy as np
import matplotlib.pyplot as plt

# --------- Green Indicators Tool --------- #
# Before using this tool make sure that the parameters used are the latest version of the parameters #

# -- Definitions -- #
class Material:
    def __init__(self, M, RC, n_rc, E_p):
        #Recycling Parameters
        self.M = M          # Part weight [kg]
        self.RC = RC        # Recyclability percentage [-]
        self.n_rc = n_rc    # Recycling energy usage / virgin production energy usage [-]
        self.E_p = E_p      # Production energy [J/kg]
        #self.n_m = n_m     # Material scarcity index (Only relevant for batteries)

def recycle(list):
    GI_1 = 0
    for i in range(0, len(list), 1):
        GIsub = (list[i].M*list[i].RC*list[i].n_rc+(1-list[i].RC)*list[i].M)*list[i].E_p
        GI_1 = GI_1 + GIsub
    print(GI_1)
    return GI_1

def Efficiency(a):
    C_D_0 = a.CD0_wing + a.CD0_tail + a.CD0_fl + a.CD0_lg
    C_D = C_D_0 + (a.CL ** 2) / (pi * a.A* a.e)
    GI_2 = a.CL / C_D * 1/(a.MTOW + a.dW * a.G) * a.n_prop * a.n_engine * a.n_pmad * a.n_cooling * a.n_fuelcell
    print("Efficiency indicator =", GI_2)
    return GI_2

def Noise(a):
    P_br = a.Pto / a.n_prop
    M_t = a.nrpm * pi * a.Dp / (60*a.c)
    GI_3 = 83.4 + 15.3 * log10(P_br/1000) - 20 * log10(a.Dp) + 38.5 * M_t - 3 * (a.B-2) + 10 * log10(a.N) - 20 * log(a.r)
    print("Noise indicator =", GI_3, "[dB]")
    return GI_3

# -- Parameter List -- #


# Efficiency and Noise inputs #
class Parameters1:
    def __init__(self):
        #General Parameters
        self.MTOW = 8618        # Maximum Take-off Weight [kg]

        #Efficiency Parameters
        self.CL = 1.3           #CL [-] (Generally use cruise conditions)
        self.CD0_lg = 0         #CD0 landing gear [-] (!!Set to 0 when not designing for takeoff!!)
        self.CD0_wing = 0.012   #CD0 wing [-]
        self.CD0_tail = 0.008   #CD0 tail [-]
        self.CD0_fl = 0.01      #CD0 fuselage [-]
        self.A = 10             #Aspect ratio [-]
        self.e = 0.8            #Oswald effiency [-]
        self.dW = -200             #Extra added weight [kg] (use negative value for weight reduction)
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

class Parameters2:
    def __init__(self):
        #General Parameters
        self.MTOW = 8618        # Maximum Take-off Weight [kg]

        #Efficiency Parameters
        self.CL = 1.3           #CL [-] (Generally use cruise conditions)
        self.CD0_lg = 0  # CD0 landing gear [-] (!!Set to 0 when not designing for takeoff!!)
        self.CD0_wing = 0.012  # CD0 wing [-]
        self.CD0_tail = 0.008  # CD0 tail [-]
        self.CD0_fl = 0.01  # CD0 fuselage [-]
        self.A = 10             #Aspect ratio [-]
        self.e = 0.8            #Oswald effiency [-]
        self.dW = 0             #Extra added weight [kg] (use negative value for weight reduction)
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

design = Parameters1()
design2 = Parameters2()

#Effiency comparison
A = Efficiency(design)
B = Efficiency(design2)

#Noise comparison
C = Noise(design)
D = Noise(design2)

# Recyclability inputs #
carbon2 = Material(100, 0.4, 0.5, 20000)
alluminium2 = Material(600, 0.8, 0.2, 1500)

carbon1 = Material(400, 0.4, 0.5, 20000)
alluminium1 = Material(100, 0.8, 0.2, 1500)

Option1 = recycle([carbon1, alluminium1])
Option2 = recycle([carbon2, alluminium2])