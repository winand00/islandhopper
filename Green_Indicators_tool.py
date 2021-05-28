from math import *
import numpy as np
import matplotlib.pyplot as plt

# --------- Green Indicators Tool --------- #
# Before using this tool make sure that the parameters used are the latest version of the parameters #

# -- Parameter List -- #

# GI1: Design efficiency #

class Concept1:
    def __init__(self):
        self.CL = 1.3 #Generally use cruise conditions
        self.C_D_0 = 0.03
        self.A = 10 #Aspect ratio
        self.e = 0.8 #Oswald
        self.MTOW = 8618 #[kg]
        self.dW = 1 #[kg] Extra added weight, use negative value for weight reduction
        self.G = 4.73 #Get this value from the Iteration_Method_tool.xcl file

def Efficiency(a):
    C_D = a.C_D_0 + (CL_value ** 2) / (pi * A_value * a.e)
    GI1 = CL/C_D * 1/(a.MTOW + a.dW * a.G) * n_prop * n_engine * n_pmad * n_cooling *

Efficiency(Concept1)