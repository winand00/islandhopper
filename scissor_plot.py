import numpy as np
from matplotlib import pyplot as plt


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
    # Set RJ85 values
    xbar_ac = 0.1038
    cbar = 3.17
    l_h = 12.35
    VhV = 1
    Cl_alpha_h = 5.9876
    Cl_alpha_Ah = 6.4334
    dedalpha = 0.3406
    Cm_ac = -0.986
    Cl_h = -0.69648
    Cl_Ah = 3.3992

    # True ShS
    ShS_true = 15.61 / 77.3

    # Cg locations RJ85 from loading diagram
    min_x_cg = 0.2908579272411666
    max_x_cg = 0.4842262602070433

    stab = stability_line(Cl_alpha_h, Cl_alpha_Ah, dedalpha, l_h, cbar, VhV, xbar_ac)
    contr = controllability_line(xbar_ac, cbar, Cl_h, Cl_Ah, l_h, VhV, Cm_ac)

    plt.plot([point[0] for point in stab], [point[1] for point in stab], label='Stability limit')
    plt.plot([point[0] for point in contr], [point[1] for point in contr], label='Controllability limit')
    plt.plot([point[0] - 0.05 for point in stab], [point[1] for point in stab], 'r--', label='Safety margin')
    plt.plot(np.linspace(min_x_cg, max_x_cg, 10), [ShS_true] * len(np.linspace(min_x_cg, max_x_cg, 10)), '--', label='RJ-85')

    plt.title('Scissor plot RJ85', weight='bold')
    plt.legend(loc='best')
    plt.grid(b=True, which='major', axis='y')
    plt.xticks(np.linspace(0, 1, 11))
    plt.yticks(np.linspace(0, 0.6, 7))
    plt.xlabel('$x_{cg}/MAC$ [-]')
    plt.ylabel('$S_{h}/S$ [-]')
    plt.show()

    # Set RJXX values
    xbar_ac = 0.1307
    Cl_alpha_Ah = 6.5949
    dedalpha = 0.30745
    Cm_ac = -0.919099
    Cl_Ah = 5.4285

    # True ShS
    ShS_true = 15.61 / 77.3

    # Cg locations RJXX from loading diagram
    min_x_cg = 0.006657003223463721
    max_x_cg = 0.3324739509679991

    stab = stability_line(Cl_alpha_h, Cl_alpha_Ah, dedalpha, l_h, cbar, VhV, xbar_ac)
    contr = controllability_line(xbar_ac, cbar, Cl_h, Cl_Ah, l_h, VhV, Cm_ac)

    plt.plot([point[0] for point in stab], [point[1] for point in stab], label='Stability limit')
    plt.plot([point[0] for point in contr], [point[1] for point in contr], label='Controllability limit')
    plt.plot([point[0] - 0.05 for point in stab], [point[1] for point in stab], 'r--', label='Safety margin')
    plt.plot(np.linspace(min_x_cg, max_x_cg, 10), [ShS_true] * len(np.linspace(min_x_cg, max_x_cg, 10)), '--',
             label='RJ-XX')

    plt.title('Scissor plot RJXX', weight='bold')
    plt.legend(loc='best')
    plt.grid(b=True, which='major', axis='y')
    plt.xticks(np.linspace(0, 1, 11))
    plt.xlabel('$x_{cg}/MAC$ [-]')
    plt.ylabel('$S_{h}/S$ [-]')
    plt.show()




