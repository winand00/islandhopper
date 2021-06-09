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
    # print(f"S_spheres: {S_spheres} \n"
    #       f"S_fillets: {S_fillets} \n"
    #       f"S_cylinders: {S_cylinders} \n"
    #       f"S_centers: {S_centers} \n"
    #       f"S_lenses: {S_lenses} \n"
    #       f"S_spheresfillets: {S_spheresfillets} \n"
    #       f"S: {S}")
    return S, S_spheres, S_fillets, S_cylinders, S_centers, S_lenses, S_spheresfillets

if __name__ == "__main__":
    result = []
    zero_point = (0, 0)

    for i in np.arange(0.2, 0.0001, -0.00001):
        result.append(multi_cell_s(2, 2, 2, i))
        if -0.0001 < multi_cell_s(2, 2, 2, i)[0] and (multi_cell_s(2, 2, 2, i)[0] < 0.0001):
            print(f"Surface area is 0 at R = {i}")
            zero_point = (i, multi_cell_s(2, 2, 2, i)[0])
    x = np.arange(0.2, 0.0001, -0.00001)
    plt.plot(x, [entry[0] for entry in result], label='$S$')
    plt.plot(x, [entry[1] for entry in result], '--', label='$S_{spheres}$')
    plt.plot(x, [entry[2] for entry in result], '--', label='$S_{fillets}$')
    plt.plot(x, [entry[3] for entry in result], '--', label='$S_{cylinders}$')
    plt.plot(x, [entry[4] for entry in result], '--', label='$S_{centers}$')
    plt.plot(x, [entry[5] for entry in result], '--', label='$S_{lenses}$')
    plt.plot(x, [entry[6] for entry in result], '--', label='$S_{spheresfillets}$')

    plt.title('Surface Area of the Hydrogen Tank', weight='bold')
    plt.legend(loc='best')
    plt.grid(b=True, which='major', axis='y')
    plt.xlabel('$Radius, R$ [m]')
    plt.ylabel('$Surface area, S$ [$m^2$]')
    plt.show()
