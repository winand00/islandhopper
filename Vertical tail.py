from math import *
import numpy as np

class Parameters:
    def __init__(self):
        self.ye = 3                 #outboard distance of engine [m]
        self.beta = 20              #sideslip [rad]
        self.dT_e = 500             #delta trust [N]
        self.CL_max = 2.8           #CLmax [-]
        self.W = 8618*9.81          #maximum take-off weight [N]
        self.lv = 7.296             #tail length [m]
        self.Sw = 45                #surface area [m^2]
        self.b = 20.12              #span [m]
        self.Peq = 750              #equivalent power [hp]
        self.Wpmax = 1900*9.81      #max payload weight [N]


def method2(a):
    Sv = 0.07*a.Sw*a.b/a.lv
    print("method2", Sv)
    return Sv

def method1(a):
    x = (a.ye/a.lv)*(a.Peq*a.CL_max)/(a.W - a.Wpmax)
    print("method1", x)
    return x
tail = Parameters()

method1(tail)
method2(tail)

