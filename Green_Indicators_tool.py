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

def recycle(a):
    GI_1 = 0
    for i in range(0, len(a), 1):
        GIsub = (a[i].M*a[i].RC*a[i].n_rc+(1-a[i].RC)*a[i].M)*a[i].E_p
        GI_1 = GI_1 + GIsub
    return GI_1

def Efficiency(a):
    S = a.S_wing + a.S_tail + a.S_fl + a.S_lg
    CL = (a.CL_wing * a.S_wing + a.CL_tail * a.S_tail + a.CL_fl * a.S_fl) / S
    C_D_0 = (a.CD0_wing * a.S_wing + a.CD0_tail * a.S_tail + a.CD0_fl * a.S_fl+ a.CD0_lg * a.S_lg) / S
    C_D_i = (a.CL_wing ** 2) / (pi * a.A_wing* a.e_wing) + (a.CL_tail ** 2) / (pi * a.A_tail* a.e_tail)
    C_D = C_D_0 + C_D_i
    GI_2 = CL / C_D * 1/(a.MTOW + a.dW * a.G) * a.n_prop * a.n_engine * a.n_pmad * a.n_cooling * a.n_fuelcell
    return GI_2

def Noise(a):
    P_br = a.Pto / a.n_prop
    M_t = a.nrpm * pi * a.Dp / (60*a.c)
    GI_3 = 83.4 + 15.3 * log10(P_br/1000) - 20 * log10(a.Dp) + 38.5 * M_t - 3 * (a.B-2) + 10 * log10(a.N) - 20 * log(a.r)
    return GI_3

def tool(a, b, c, d):
    print("")
    print("----- Comparison between design options -----")
    print("G1 Recyclability indicator:")
    print("Option 1 = ", recycle(b), "[J]")
    print("Option 2 = ", recycle(d), "[J]")
    print("Option 1 is ", recycle(b)/recycle(d), "relative to Option 2")
    print("")
    print("G2 Efficiency indicator:")
    print("Option 1 = ", Efficiency(a))
    print("Option 2 = ", Efficiency(c))
    print("Option 1 is ", Efficiency(a) / Efficiency(c), "relative to Option 2")
    print("")
    print("G3 Noise indicator: ")
    print("Option 1 = ", Noise(a), "[dB]")
    print("Option 2 = ", Noise(c), "[dB]")
    print("Option 1 is ", Noise(a) / Noise(c), "relative to Option 2")


# -- Parameter List -- #

# Efficiency and Noise inputs #
class Parameters1:
    def __init__(self):
        #Efficiency Indicator Parameters option 1
        self.CL_wing = 1.3              # CL  wing [-] (Generally use cruise conditions)
        self.CD0_wing = 0.012           # CD0 wing [-]
        self.S_wing = 100               # Surface area wing [m^2]
        self.A_wing = 10                # Aspect ratio [-]
        self.e_wing = 0.8               # Oswald effiency wing [-]

        #Tail option 1
        self.CL_tail = 1.3              # CL  tail [-] (Generally use cruise conditions)
        self.CD0_tail = 0.008           # CD0 tail [-]
        self.S_tail = 50                # Surface area tail [m^2]
        self.A_tail = 10                # Aspect ratio [-]
        self.e_tail = 0.8               # Oswald effiency tail [-]

        #Fuselage option 1
        self.CL_fl = 1.3                # CL fuselage [-] (Generally use cruise conditions)
        self.CD0_fl = 0.01              # CD0 fuselage [-]
        self.S_fl = 10                  # Surface area fuselage  [m^2]

        #Landing Gear option 1
        self.CD0_lg = 0                 #CD0 landing gear [-] (!!Set to 0 when not designing for takeoff!!)
        self.S_lg = 1                   #Surface area landing gear  [m^2]

        # General Parameters option 1
        self.MTOW = 8618                # Maximum Take-off Weight [kg]
        self.dW = -200                  #Extra added weight [kg] (use negative value for weight reduction)
        self.G = 4.73                   #Weight Growth Factor [-] (Get this value from the Iteration_Method_tool.xcl file)
        self.n_prop = 0.95              #Propeller efficiency [-]
        self.n_engine = 0.93            #Engine efficiency [-]
        self.n_pmad = 0.9               #PMAD efficiency [-]
        self.n_cooling = 0.4            #Cooling efficiency [-]
        self.n_fuelcell = 0.65          #Fuel cell efficiency [-]

        #Noise Parameters option 1
        self.Pto = 1100000              #Take-off Power [W]
        self.Dp = 3.0                   #Propeller diameter [m]
        self.B = 5                      #Number of blades per propeller [-]
        self.nrpm = 2500                #Engine RPM [evolutions/minute]
        self.c = 343                    #Speed of sound [m/s]
        self.N = 2                      #Number of propellers [-]
        self.r = 350                    #Distance to observer [m]

class Parameters2:
    def __init__(self):
        # Efficiency Indicator Parameters option 2
        self.CL_wing = 1.3              # CL  wing [-] (Generally use cruise conditions)
        self.CD0_wing = 0.012           # CD0 wing [-]
        self.S_wing = 100               # Surface area wing [m^2]
        self.A_wing = 10                # Aspect ratio [-]
        self.e_wing = 0.8               # Oswald effiency wing [-]

        # Tail option 2
        self.CL_tail = 1.3              # CL  tail [-] (Generally use cruise conditions)
        self.CD0_tail = 0.008           # CD0 tail [-]
        self.S_tail = 50                # Surface area tail [m^2]
        self.A_tail = 10                # Aspect ratio [-]
        self.e_tail = 0.8               # Oswald effiency tail [-]

        # Fuselage option 2
        self.CL_fl = 1.3                # CL fuselage [-] (Generally use cruise conditions)
        self.CD0_fl = 0.01              # CD0 fuselage [-]
        self.S_fl = 10                  # Surface area fuselage  [m^2]

        # Landing Gear option 2
        self.CD0_lg = 0                 # CD0 landing gear [-] (!!Set to 0 when not designing for takeoff!!)
        self.S_lg = 1                   # Surface area landing gear  [m^2]

        # General Parameters option 2
        self.MTOW = 8618                # Maximum Take-off Weight [kg]
        self.dW = -200                  # Extra added weight [kg] (use negative value for weight reduction)
        self.G = 4.73                   # Weight Growth Factor [-] (Get this value from the Iteration_Method_tool.xcl file)
        self.n_prop = 0.95              # Propeller efficiency [-]
        self.n_engine = 0.93            # Engine efficiency [-]
        self.n_pmad = 0.9               # PMAD efficiency [-]
        self.n_cooling = 0.4            # Cooling efficiency [-]
        self.n_fuelcell = 0.65          # Fuel cell efficiency [-]

        # Noise Parameters option 2
        self.Pto = 1100000              # Take-off Power [W]
        self.Dp = 3.0                   # Propeller diameter [m]
        self.B = 5                      # Number of blades per propeller [-]
        self.nrpm = 2500                # Engine RPM [evolutions/minute]
        self.c = 343                    # Speed of sound [m/s]
        self.N = 2                      # Number of propellers [-]
        self.r = 350                    # Distance to observer [m]

design1 = Parameters1()
design2 = Parameters2()


# Recyclability inputs #
carbon2 = Material(100, 0.4, 0.5, 20000)
alluminium2 = Material(600, 0.8, 0.2, 1500)

carbon1 = Material(400, 0.4, 0.5, 20000)
alluminium1 = Material(100, 0.8, 0.2, 1500)

Option1 = [carbon1, alluminium1]
Option2 = [carbon2, alluminium2]


tool(design1, Option1, design2, Option2)
