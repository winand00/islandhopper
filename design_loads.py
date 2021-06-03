from math import *
import numpy as np
from matplotlib import pyplot as plt

def design_loads():
    # Aircraft parameters
    S = 45  # [m2]
    MAC = 2.24   # [m]
    Cl_alpha = 4    # ??? Cl_alpha [rad-1]
    CL_clean = 2   # ???
    CL_flaps = 2.4  # ???

    g = 9.80665  # [m/s2]
    W = 8618.25503 * 9.80665  # [N]
    rho_zero = 1.225  # sealevel density [kg/m3]
    rho = 0.6527  # [kg/m3] density from 6096-0 [m]

    VC = 31.9426 * sqrt(19000 / (S*10.76391)) * 0.514444 # W/S in lbs/ft2, Value of 31 changes with W/S, see CS-23
    VD = 1.3880 * VC            # Value of 1.38 changes with W/S, see CS-23
    Ude_C = 15.24  # [m/s]
    Ude_D = 7.62  # [m/s]
    Ude_B = 20.12  # [m/s]

    # Aeroplane mass ratio formula
    mug = (2*W/S)/(rho*MAC*Cl_alpha*g)
    # Gust allevation factor
    kg = (0.88*mug)/(5.3 + mug)

    # Gust load factor at VC
    delta_n_C = (kg*rho_zero*Ude_C*VC*Cl_alpha)/(2*W/S)
    n_pos_C = 1 + delta_n_C
    n_neg_C = 1 - delta_n_C
    """
    print('At VC:')
    print('n_pos = ', n_pos_C)
    print('n_neg = ', n_neg_C)
    print(' ')
    """
    # Gust load factor at VD
    delta_n_D = (kg * rho_zero * Ude_D * VD * Cl_alpha) / (2 * W / S)
    n_pos_D = 1 + delta_n_D
    n_neg_D = 1 - delta_n_D
    """
    print('At VD:')
    print('n_pos = ', n_pos_D)
    print('n_neg = ', n_neg_D)
    print(' ')
    """
    # Calculation of VB, first using VS1 sqrt ng
    VS1 = sqrt(W / (0.5 * rho * S * CL_clean))  # [m/s]
    VB_1 = sqrt(n_pos_C) * VS1

    # Second way of VB calculation, using intersection between stall curve and VB air gust line
    V = np.arange(0, 130, 1)
    n_B_line = 1 + (kg * rho_zero * Ude_B * V * Cl_alpha) / (2 * W / S)
    n_S_line = 0.5 * rho * V ** 2 * CL_clean / (W/S)  # Assumed CLmax clean is 2
    for i in np.arange(0, 130, 1):
        diff = n_B_line[i] - n_S_line[i]
        if diff > 0 and diff < 0.03:
            VB_2 = V[i]

    # Choose lowest VB, cannot be higher than VC
    if VB_1 < VB_2:
        VB = VB_1
    else:
        VB = VB_2
    if VB > VC:
        VB = VC


    # Gust load factor at VB
    delta_n_B = (kg * rho_zero * Ude_B * VB * Cl_alpha) / (2 * W / S)
    n_pos_B = 1 + delta_n_B
    n_neg_B = 1 - delta_n_B
    """
    print('At VB:')
    print('n_pos = ', n_pos_B)
    print('n_neg = ', n_neg_B)
    """

    # Gust load factor with flaps out
    VS2 = sqrt(W / (0.5 * rho * S * CL_flaps))  # [m/s]
    VF_1 = 1.4 * VS1
    VF_2 = 1.8 * VS2
    if VF_1 < VF_2:
        VF = VF_1
    else:
        VF = VF_2
    delta_n_F = (kg * rho_zero * Ude_D * VF * Cl_alpha) / (2 * W / S)
    n_pos_F = 1 + delta_n_F

    # Calculate max and minimum load factors
    lf_pos = max(n_pos_C, n_pos_D, n_pos_B, 2.9278)
    lf_neg = min(n_neg_C, n_neg_D, n_neg_B, -1.171)
    lf_pos_flaps = max(n_pos_F, 2)
    return lf_pos, lf_neg, lf_pos_flaps


if __name__ == "__main__":
    lf_pos, lf_min = design_loads()
    print(lf_pos, lf_min)







