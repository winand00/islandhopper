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

    # Tail parameters
    aht = 2  #???    # Tail lift curve slope
    Sht = 5   #???   # Tail surface area
    deda = 0.7   #???    # Lift downwash per alpha

    g = 9.80665  # [m/s2]
    W = 8618.25503 * 9.80665  # [N]
    rho_zero = 1.225  # sealevel density [kg/m3]
    rho = 0.7708#0.6527  # [kg/m3] density from 6096-0 [m]

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
    #delta_n_C_check = kg*Ude_C*3.28*VC*1.944*Cl_alpha/(498*(0.020885*W/S))
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
    n_S_line = 0.5 * rho * V ** 2 * CL_clean / (W/S)
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

    # Change in tail load due to gusts
    #delta_L_F = kg * Ude_D*3.2808 * VF*1.9439 * aht * Sht*10.7639 / 498 * (1 - deda) * 4.44822
    #delta_L_C = kg * Ude_C*3.2808 * VC*1.9439 * aht * Sht*10.7639 / 498 * (1 - deda) * 4.44822
    #delta_L_D = kg * Ude_D*3.2808 * VD*1.9439 * aht * Sht*10.7639 / 498 * (1 - deda) * 4.44822

    #delta_L = max(delta_L_F, delta_L_C, delta_L_D)

    # Calculate max and minimum load factors
    lf_pos = max(n_pos_C, n_pos_D, n_pos_B, 2.9278)
    lf_neg = min(n_neg_C, n_neg_D, n_neg_B, -1.171)
    lf_pos_flaps = max(n_pos_F, 2)

    # Plot loading diagram
    fig, ax = plt.subplots()

    # Get 2% margin lines:

    # Manoeuvre lines
    VA = sqrt(2.9278 * W / (S * rho * 0.5 * CL_clean))
    V_VA = np.arange(0, VA + 0.6, 1)
    n_VA = CL_clean / (W/S) * 0.5 * rho * V_VA **2
    plt.plot(V_VA,n_VA,color='k')
    V_VD = np.arange(VA, VD, 1)
    n_VD = 2.9278 * np.ones(len(V_VD))
    plt.plot(V_VD,n_VD,color='k')
    V_VS = np.arange(0, VS1 + 1, 1)
    n_VS = - CL_clean / (W/S) * 0.5 * rho * V_VS **2 * 1.16
    plt.plot(V_VS,n_VS,color='k')
    V_VC = np.arange(VS1, VC, 1)
    n_VC = -1.171 * np.ones(len(V_VC))
    plt.plot(V_VC,n_VC,color='k')
    V_VD = np.arange(VC, VD + 1, 1)
    n_VD = -1.171 + 1.171 / (VD-VC) * (V_VD - VC)
    plt.plot(V_VD,n_VD,color='k')

    # Gust lines
    V_VB = np.arange(0, VB, 1)
    n_VB = 1 + (n_pos_B -1) / VB * V_VB
    plt.plot(V_VB, n_VB, color='r')
    V_VBVC = np.arange(VB, VC, 1)
    n_VBVC = n_pos_B + (n_pos_C - n_pos_B) / (VC - VB) * (V_VBVC - VB)
    plt.plot(V_VBVC, n_VBVC, color='r')
    n_VCVD = n_pos_C + (n_pos_D - n_pos_C) / (VD - VC) * (V_VD - VC)
    plt.plot(V_VD, n_VCVD, color='r')
    n_VB_neg = 1 + (n_neg_B - 1) / VB * V_VB
    plt.plot(V_VB, n_VB_neg, color='r')
    n_VBVC_neg = n_neg_B + (n_neg_C - n_neg_B) / (VC - VB) * (V_VBVC - VB)
    plt.plot(V_VBVC, n_VBVC_neg, color='r')
    n_VCVD_neg = n_neg_C + (n_neg_D - n_neg_C) / (VD - VC) * (V_VD - VC)
    plt.plot(V_VD, n_VCVD_neg, color='r')
    plt.vlines(x=VD, ymin=0, ymax=2.9278, colors='k')
    #plt.axvline(VD, ymin=0, ymax=2.9278, color='black', linestyle='-')
    plt.title('Load Diagram Hopper')
    plt.xlabel('V [m/s]')
    plt.ylabel('n [-]')
    plt.ylim((-1.5, 3.2))
    plt.xlim(0, VD + 10)
    plt.xticks(np.arange(0, VD + 10, 10))
    plt.yticks(np.arange(-1.5, 3.2, 0.5))
    ax.set(facecolor='w')
    # Legend:
    #plt.legend(loc=2, prop={'size': 8})
    plt.grid(axis='y')

    return lf_pos, lf_neg, lf_pos_flaps


def tail_load_elevator():
    delta_n = lf_pos             # load factor increment
    M = 8618.25503  # [kg]
    g = 9.80665
    x_cg = 6
    x_ac = 5.5
    x_cg_ac = x_cg - x_ac # ???   # [m] distance from ac to cg
    l_t = 9  # ???         # tail arm
    S_h_t = 6   # ???      # hor tail area
    S = 45                 # wing surface
    aht = 2   # ???        # lift curve slope hor tail
    a = 3     # ???        # lift curve slope wing
    deda = 0.7             # downwash change with alpha
    rho_zero = 1.225       # density sealevel

    # Change in tail load due to elevator deflection
    delta_P = delta_n * M * g * ((x_cg_ac/l_t) - (S_h_t/S) * ((aht/a) * (1-deda) - rho_zero/2 * (S_h_t * aht * l_t / M)))

    return delta_P

if __name__ == "__main__":
    lf_pos, lf_min, v = design_loads()
    fsadf = tail_load_elevator()
    plt.show()







