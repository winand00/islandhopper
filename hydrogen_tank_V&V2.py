import numpy as np
import math
from matplotlib import pyplot as plt

# CONSTANTS
k_overlap = 0.75
P = 5.7

t_ply = 0.00014

t_metal_liner = 0.001
rho_metal_liner = 2700

t_polyamide_liner = 0.001
rho_polyamide_liner = 1010

rho_comp = 1800
k_comp = 60

rho_ins = 32
k_ins = 0.022

y = 1.4
R_over_t = 130
s_over_R = 0.35           #(fig 4.22)
r = 0.1                   #(fig 4.23a) [R_fillet/R]
t_junction_over_R = 0.027 #(fig 4.23d)
phi = 24                  #degrees (fig 4.23b) [layup angle]
i_ratio = 2.5             #(fig 4.23b)   [ratio between 0 degrees fibers and phi degrees fibres]
k_overlap = 0.75          #(fig 5.5)

rho_hydr = 70 # [kg/m3]
d_H_vap = 446100 #J/kg
BOR_percentage = 0.016 # [-]
d_T = 293.15 # temp diff between 20 K and 40 C

def multi_cell_v(m, n, p, R):
    """"Function to calculate multi cell volume in m3."""
    d = y * R
    R_fillet = r * R
    s = s_over_R * R
    t_junction = t_junction_over_R * R

    N_cells = m * n * p
    V_spheres = N_cells * 4 / 3 * np.pi * R ** 3

    N_junctions = (n * (m - 1) + m * (n - 1)) * p + m * n * (p - 1)
    N_centers = ((m - 1) * (n - 1)) * p + (p - 1) * (n - 1)

    V_lensjunction = (np.pi * (2 * R - d) ** 2 * (d ** 2 + 4 * d * R)) / (12 * d)
    V_lenses = N_junctions * V_lensjunction

    # R_fillet = R*r

    theta_1 = np.arcsin((d / 2) / (R + R_fillet))  # e1 4.42
    h_inter = R_fillet * np.sin(theta_1)
    R_ring = np.sqrt(R ** 2 - (d / 2) ** 2)  # eq 4.41

    R_ringfinal = R_ring + ((R * np.cos(theta_1) - R_ring) - (h_inter * np.tan(theta_1 / 2)))
    V_centerlens = (np.pi * (2 * R_ringfinal - d) ** 2 * (d ** 2 + 4 * d * R_ringfinal)) / (12 * d)
    V_centers = N_centers * V_centerlens

    h_center = (2 * R - d * 2 ** 0.5) / 2
    r_cyl = max((4 * s * k_overlap) / 2 / np.pi, 2 * t_junction, h_center)
    l_cyl = np.sqrt(R_ringfinal ** 2 - (d / 2 - r_cyl) ** 2)
    V_cylinder = np.pi * r_cyl ** 2 * l_cyl
    V_cylinders = N_centers * V_cylinder

    A_triangle = ((2 * R_fillet * np.sin(theta_1)) * (R * np.cos(theta_1) - R_ring)) / 2
    A_cap = R_fillet ** 2 / 2 * (2 * theta_1 - np.sin(2 * theta_1))
    V_fillet = (A_triangle - A_cap) * (2 * np.pi * R_ring)
    V_fillets = N_junctions * V_fillet

    V = (V_spheres + V_centers + V_fillets - V_lenses - V_cylinders) #* 1000
    return V

def multi_cell_s(m, n, p, R):
    """"Function to calculate multi cell surface area in m2."""
    d = y * R
    R_fillet = r * R
    s = s_over_R * R
    t_junction = t_junction_over_R * R

    N_cells = m * n * p
    theta_1 = np.arcsin((d / 2) / (R + R_fillet))  # e1 4.42
    h_inter = R_fillet * np.sin(theta_1)
    R_ring = np.sqrt(R ** 2 - (d / 2) ** 2)  # eq 4.41
    R_ringfinal = R_ring + ((R * np.cos(theta_1) - R_ring) - (h_inter * np.tan(theta_1 / 2)))
    N_junctions = (n * (m - 1) + m * (n - 1)) * p + m * n * (p - 1)
    N_centers = ((m - 1) * (n - 1)) * p + (p - 1) * (n - 1)
    h_center = (2 * R - d * 2 ** 0.5) / 2
    r_cyl = max((4 * s * k_overlap) / 2 / np.pi, 2 * t_junction, h_center)
    l_cyl = np.sqrt(R_ringfinal ** 2 - (d / 2 - r_cyl) ** 2)

    S_spheres = N_cells * 4 * np.pi * R ** 2
    h_cap = R - np.sqrt(R ** 2 - R_ring ** 2)
    S_lensjunction = 2 * (2 * h_cap)
    S_lenses = N_junctions * S_lensjunction

    gamma_cyl = np.arcsin(l_cyl / (2 * R_ringfinal))
    S_torus = 2 * np.pi * s * R_ring * ((2 * np.pi - gamma_cyl) / (2 * np.pi))
    S_fillets = N_junctions * S_torus
    S_spherefillet = 2 * (2 * np.pi * R * h_inter)
    S_spheresfillets = N_junctions * S_spherefillet

    S_center = 2 * (2 * np.pi * R * h_center)
    S_centers = N_centers * S_center

    S_cylinder = 2 * np.pi * (r_cyl * l_cyl + r_cyl ** 2)
    S_cylinders = N_centers * S_cylinder

    S = S_spheres + S_fillets + S_cylinders + S_centers - S_lenses - S_spheresfillets

    t = max(R / 130, 6 * t_ply)
    M_structural = rho_comp * t * (S_spheres + S_centers - S_lenses - S_spheresfillets) + rho_comp * t_junction * \
                   (S_fillets + S_cylinders)
    return S, M_structural, t


def multi_cell_m_incl_insulation(m, n, p, R, Mass):
    """"Function to calculate multi cell mass in kg."""
    S, M_structural, t = multi_cell_s(m, n, p, R)

    volume = multi_cell_v(m, n, p, R)
    tanks_needed = math.ceil(Mass/(volume * rho_hydr))

    optimum = False
    S_cryo = S

    threshold = 0.001
    t_result = []
    while not optimum:

        S_old = S_cryo

        Q = (Mass * BOR_percentage) * d_H_vap
        U = Q / (S_cryo * d_T)

        h_out = Q/3600 /S / (-1*d_T) #-4.5 #


        # t_ins = (1 / U - t / k_comp - 1 / h_out) * k_ins
        # t_result.append(t_ins)
        # print(t_ins)

        t_ins = (1 / U - t / k_comp - 1 / h_out) * k_ins
        if math.isnan(t_ins):
            break

        t_result.append(t_ins)
        print(t_ins)

        # except FloatingPointError:
        #     print('Runtime warning')
        #     break

        # try:
        #     S_cryo, M_throwaway, t_throwaway = multi_cell_s(m, n, p, R + t_ins)
        # except RuntimeWarning:
        #     print('Runtime warning')
        #     break

        S_cryo, M_throwaway, t_throwaway = multi_cell_s(m, n, p, R + t_ins)

        delta_t = abs(S_old - S_cryo)

        if delta_t < threshold:
            optimum = True

    # optimum = False
    # t_ins = 1
    #
    # threshold = 0.00001
    # while not optimum:
    #
    #     print("t_ins", t_ins)
    #
    #
    #     t_old = t_ins
    #     S_cryo, M_throwaway, t_throwaway = multi_cell_s(m, n, p, R + t_ins)
    #     Q = (Mass * BOR_percentage)  * d_H_vap
    #     U = Q / (S_cryo * d_T)
    #     t_ins = (1 / U - t / k_comp - 1 / h_out) * k_ins
    #     delta_t = abs(t_old - t_ins)
    #
    #     print("delta", delta_t)
    #
    #
    #     if delta_t < threshold:
    #         optimum = True


    M_insulation = rho_ins * t_ins * S_cryo
    M_total = M_structural + M_insulation
    return M_total, M_structural, M_insulation, S, S_cryo, t_ins, tanks_needed, t_result

if __name__ == '__main__':
    result = []

    for i in np.arange(100, 10, -1):
        print(f'START with i = {i}')
        result.append(multi_cell_m_incl_insulation(2, 2, 2, 1, i))
        print('STOP')

    # for i in range(80, 89):
    #     plt.plot(np.arange(0, len(result[i][-1])), result[i][-1])

    plt.figure()
    plt.plot(np.arange(0, len(result[80][-1])), result[80][-1], label='20 kg')
    plt.plot(np.arange(0, len(result[81][-1])), result[81][-1], label='19 kg')
    plt.plot(np.arange(0, len(result[82][-1])), result[82][-1], label='18 kg')
    plt.plot(np.arange(0, len(result[83][-1])), result[83][-1], label='17 kg')
    plt.plot(np.arange(0, len(result[84][-1])), result[84][-1], label='16 kg')
    plt.plot(np.arange(0, len(result[85][-1])), result[85][-1], label='15 kg')
    plt.plot(np.arange(0, len(result[86][-1])), result[86][-1], label='14 kg')
    plt.plot(np.arange(0, len(result[87][-1])), result[87][-1], label='13 kg')
    plt.plot(np.arange(0, len(result[88][-1])), result[88][-1], label='12 kg')
    # plt.plot(np.arange(0, len(result[89][-1])), result[89][-1], label='11 kg')

    plt.title('Insulation Thickness after Iteration (Conversion)', weight='bold')
    plt.legend(loc='best', title='m=n=p=2\nR=1\nHydrogen mass')
    plt.grid(b=True, which='major', axis='y')
    plt.xlabel('Iterations [$-$]')
    plt.ylabel('Insulation thickness, $t_{ins}$ [$m$]')
    plt.show()

    plt.figure()
    plt.plot(np.arange(0, len(result[89][-1])), result[89][-1], label='11 kg')
    plt.title('Insulation Thickness after Iteration (Diversion)', weight='bold')
    plt.legend(loc='best', title='m=n=p=2\nR=1\nHydrogen mass')
    plt.grid(b=True, which='major', axis='y')
    plt.xlabel('Iterations [$-$]')
    plt.ylabel('Insulation thickness, $t_{ins}$ [$m$]')
    plt.yscale('log')
    plt.show()










    # print(multi_cell_m_incl_insulation(2, 2, 2, 1, 12)) # (2, 2, 2, 1, 12) lowest h possible