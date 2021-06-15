import numpy as np
from matplotlib import pyplot as plt
from CGloadingdiagram import cgboundaries
from Class2_weight_estimation import cg_wing, cg_horizontal_tail, L_fuselage_whole, S_w, CL_max, V_stall, B_w


#Define functions
def stability_line(Cl_alpha_h, Cl_alpha_Ah, dedalpha, l_h, cbar, VhV, xbar_ac):
    """Returns stable S_h/S as a function of x̄_cg (location of the center of gravity divided by c̄, the mean chord)."""
    result = []
    xbar_cg = np.linspace(0, 1, 100)
    for i in xbar_cg:
        ShS = (1 / ((Cl_alpha_h / Cl_alpha_Ah) * (1 - dedalpha) * (l_h / cbar) * (VhV ** 2))) * i - (xbar_ac - 0.05) / \
              ((Cl_alpha_h / Cl_alpha_Ah) * (1 - dedalpha) * (l_h / cbar) * (VhV ** 2))
        result.append((i, ShS))
    return [point for point in result if point[1] >= 0]


def controllability_line(xbar_ac, cbar, Cl_h, Cl_Ah, l_h, VhV, Cm_ac):
    """Returns controllable S_h/S as a function of x̄_cg (location of the center of gravity divided by c̄, the mean chord)."""
    result = []
    xbar_cg = np.linspace(0, 1, 100)
    for i in xbar_cg:
        ShS = (1 / ((Cl_h / Cl_Ah) * (l_h / cbar) * (VhV ** 2))) * i + ((Cm_ac / Cl_Ah) - xbar_ac) / ((Cl_h / Cl_Ah) *
                                                                                                      (l_h / cbar) * (VhV ** 2))
        result.append((i, ShS))
    return [point for point in result if point[1] >= 0]

def get_ac(cg_wing,cbar,Cl_alpha_Ah,bf,hf,S,mac,lambda_taper,lambda_a):
    lfn = cg_wing - 0.5 * cbar
    ac_fus1 = (-1.8 / Cl_alpha_Ah) * (bf * hf * lfn) / (S * mac)
    ac_fus2 = (0.273 / (1 + lambda_taper)) * (bf * (S / b) * (b - bf)) / ((S / b) ** 2 * (b + 2.15 * bf)) * lambda_a
    xbar_ac = 0.25 * mac + ac_fus1 + ac_fus2  # (5/12.5)*lf/mac #aerodynamic centre location
    return xbar_ac


if __name__ == '__main__':
    # Set Hopper values, iteration inputs
    bf = 2.11
    hf =   2.4     #Fuselage height
    cbar = 2

#------------------Not for iteration---------------
    b = B_w
    V_low = 1.1* V_stall
    S = S_w
    lf = L_fuselage_whole  # fuselage length
    S_net = S - bf * 3  # S minus b_f * chord at root
    mac = 2.074             #mac = 2/3*0.5*cbar*((1+0.5+0.5**2)/(1+0.5))
    VhV = 0.95
    Ah = 6.73
    A = 10
    Cm0_airfoil = -0.0232
    CL0_airfoil = .2885
    eta = 0.95
    lambda_taper = 0.5
    lambda_h = 0
    lambda_a = 0

    #------------Calculated values------------
    M = 0.275
    Mh = 0.275 * VhV
    Mlow = V_low /  327.27
    beta = np.sqrt(1 - M * M)
    betah = np.sqrt(1 - Mh * Mh)
    betalow = np.sqrt(1 - Mlow * Mlow)
    SnS = S_net / S

    Cl_alpha_w = 2*np.pi*A/(2+np.sqrt(4+(A*beta/eta)**2*(1+np.tan(lambda_a)**2/beta**2)))
    print("Cl_alpha_w = ", Cl_alpha_w)
    Cl_alpha_h = 2*np.pi*Ah/(2+np.sqrt(4+(Ah*betah/eta)**2*(1+np.tan(lambda_h)**2/betah**2)))
    Cl_alpha_Ah =  Cl_alpha_w*(1+2.15*bf/b)*SnS +np.pi/2*bf**2 /S
    print("Cl_alpha_h = ", Cl_alpha_h)
    print("Cl_alpha_Ah = ", Cl_alpha_Ah)
    Cl_alpha_A = Cl_alpha_Ah + Cl_alpha_h
    print("Cl_alpha_A = ", Cl_alpha_A)

    xbar_ac = get_ac(cg_wing, cbar, Cl_alpha_Ah, bf, hf, S, mac, lambda_taper, lambda_a)
    print("x_ac=",cg_wing - cbar/2 + xbar_ac)

    l_h = cg_horizontal_tail - (cg_wing - mac/2 + xbar_ac)  # 5.93 #(5.7/12.5)*lf
    r = l_h/(b/2)
    mtv = 0.2
    Kel = 0.1124/r**2 + 0.1024/r + 2
    Kel0 = 0.1124/r**2 + 0.1024/r + 2
    T1 = (r/(r**2+mtv**2))*0.4876/(np.sqrt(r**2+0.6319+mtv**2))
    T2 = (1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2)))

    
    dedalpha = Kel/Kel0*(T1+T2)*Cl_alpha_w/3.1415926535/A   #Validated with slide formula 4/(A+2)
    print("De/dalpha is: ", dedalpha)

    #Cm_ac calculation
    Cl_alpha_w_low = 2 * np.pi * A / (2 + np.sqrt(4 + (A * betalow / eta) ** 2 * (1 + np.tan(lambda_a) ** 2 / betalow ** 2)))
    Cl_alpha_Ah_low = Cl_alpha_w_low*(1+2.15*bf/b)*SnS +np.pi/2*bf**2 /S
    Cm_acw = Cm0_airfoil * (A*np.cos(lambda_a)*np.cos(lambda_a))/(A+2*np.cos(lambda_a))
    dfC_m_ac = - 0.33       #From sead lecture 5 slide 21
    dfusCm_ac = -1.8 * (1-2.5*bf/lf) * (np.pi * bf * hf * lf)/(4*S*cbar) * (CL0_airfoil / Cl_alpha_Ah_low)
    Cm_ac = Cm_acw +  dfC_m_ac  + dfusCm_ac        #-0.986
    print("Cm_ac=", Cm_ac)
    Cl_h = -0.35*Ah**0.333
    Cl_Ah = CL_max - Cl_h

    # True ShS
    ShS_true = 7 / 40

    # Cg locations Hopper from loading diagram
    min_x_cg, max_x_cg = cgboundaries() 

    stab = stability_line(Cl_alpha_h, Cl_alpha_Ah, dedalpha, l_h, cbar, VhV, xbar_ac)
    contr = controllability_line(xbar_ac, cbar, Cl_h, Cl_Ah, l_h, VhV, Cm_ac)

    plt.plot([point[0] for point in stab], [point[1] for point in stab], label='Stability limit')
    plt.plot([point[0] for point in contr], [point[1] for point in contr], label='Controllability limit')
    plt.plot([point[0] - 0.05 for point in stab], [point[1] for point in stab], 'r--', label='Safety margin')
    plt.plot(np.linspace(min_x_cg, max_x_cg, 10), [ShS_true] * len(np.linspace(min_x_cg, max_x_cg, 10)), '--', label='Hopper')

    plt.title('Scissor plot Hopper', weight='bold')
    plt.legend(loc='best')
    plt.grid(b=True, which='major', axis='y')
    plt.xticks(np.linspace(0, 1, 11))
    plt.yticks(np.linspace(0, 0.9, 10))
    plt.xlabel('$x_{cg}/MAC$ [-]')
    plt.ylabel('$S_{h}/S$ [-]')
    plt.show()
