from math import *
import numpy as np

def design_loads():
    # Constant parameters
    W = 8618.25503*9.80665 # [N]
    rho_zero = 1.225  # sealevel density [kg/m3]
    S = 45  # [m2]
    MAC = 2.24   # [m]
    g = 9.80665  # [m/s2]
    Cl_alpha = 4    # ??? Cl_alpha [rad-1]

    # Variable parameters
    rho = 1.225#0.6527 # [kg/m3] density from 6096-0 [m]

    VC = 31.9426 * sqrt(19000 / (S*10.76391)) * 0.514444 # W/S in lbs/ft2, Value of 31 changes with W/S, see CS-23
    VD = 1.3880 * VC   # Value of 1.38 changes with W/S, see CS-23
    Ude_C = 15.24  # [m/s]
    Ude_D = 7.62  # [m/s]
    Ude_B = 20.12  # [m/s]

    # Aeroplane mass ratio formula
    mug = (2*W/S)/(rho*MAC*Cl_alpha*g)
    # Gust allevation factor
    kg = (0.88*mug)/(5.3 + mug)
    # Gust load factor at VC
    delta_n = (kg*rho_zero*Ude_C*VC*Cl_alpha)/(2*W/S)
    n_pos_C = 1 + delta_n
    n_neg_C = 1 - delta_n
    print('At VC:')
    print('n_pos = ', n_pos_C)
    print('n_neg = ', n_neg_C)
    print(' ')

    # Gust load factor at VC
    delta_n = (kg * rho_zero * Ude_D * VD * Cl_alpha) / (2 * W / S)
    n_pos_D = 1 + delta_n
    n_neg_D = 1 - delta_n
    print('At VD:')
    print('n_pos = ', n_pos_D)
    print('n_neg = ', n_neg_D)

    print(VC)
    print(VD)


if __name__ == "__main__":
    design_loads()







