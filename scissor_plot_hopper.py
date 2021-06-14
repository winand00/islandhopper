import numpy as np
from matplotlib import pyplot as plt
from CGloadingdiagram import cgboundaries


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


if __name__ == '__main__':
    # Set Hopper values
    M = 0.274
    beta = np.sqrt(1-M*M)
    eta = 0.95
    lambda_h = 0
    lambda_a = 0
    b = 20.12
    bf = 2.09
    lf = 12.24 #fuselage length
    S = 45
    SnS = 0.9
    cbar = 2.24
    mac = 2/3*0.5*cbar*((1+0.5+0.5**2)/(1+0.5))
    print("Mac:",mac)
    xbar_ac = 0.25*mac #(5/12.5)*lf/mac #aerodynamic centre location
    l_h = (5.7/12.5)*lf
    VhV = 0.95
    Ah = 6.73
    A = 9
    
    Cl_alpha_w = 2*np.pi*A/(2+np.sqrt(4+(A*beta/eta)**2*(1+np.tan(lambda_a)**2/beta**2)))
    print("Cl_alpha_w = ", Cl_alpha_w)
    Cl_alpha_h = 2*np.pi*Ah/(2+np.sqrt(4+(Ah*beta/eta)**2*(1+np.tan(lambda_h)**2/beta**2)))
    Cl_alpha_Ah =  5.8 #Cl_alpha_w*(1+2.15*bf/b)*SnS +np.pi/2*bf**2/S
    print("Cl_alpha_h = ", Cl_alpha_h)
    print("Cl_alpha_Ah = ", Cl_alpha_Ah)
    Cl_alpha_A = Cl_alpha_Ah + Cl_alpha_h
    print("Cl_alpha_A = ", Cl_alpha_A)
    
    r = 0.867
    mtv = 0.067
    Kel = 1
    Kel0 = 1
    T1 = (r/(r**2+mtv**2))*0.4876/(np.sqrt(r**2+0.6319+mtv**2))
    T2 = (1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2)))
    
    
    dedalpha = Kel/Kel0*(T1+T2)*Cl_alpha_w/3.1415926535/A
    Cm_ac = -0.986 #!
    Cl_h = -0.5*Ah**0.333
    Cl_Ah = 2.4- Cl_h #!

    # True ShS
    ShS_true = 13 / 45

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
    plt.yticks(np.linspace(0, 0.6, 7))
    plt.xlabel('$x_{cg}/MAC$ [-]')
    plt.ylabel('$S_{h}/S$ [-]')
    plt.show()
