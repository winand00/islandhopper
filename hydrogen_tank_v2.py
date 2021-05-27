import numpy as np
import math

# CONSTANTS
k_overlap = 0.75      #Used in thesis
P = 5.7
# sigma_allowable = 1
t_ply = 0.00014       #[m] #Table 4-3
rho_comp = 2700
k_comp = 60
rho_ins = 200
k_ins = 0.025
y = 1.4
R_over_t = 130
s_over_R = 0.35           #(fig 4.22)
r = 0.1                   #(fig 4.23a) [R_fillet/R]
t_junction_over_R = 0.027 #(fig 4.23d)
phi = 24                  #degrees (fig 4.23b) [layup angle]
i_ratio = 2.5             #(fig 4.23b)   [ratio between 0 degrees fibers and phi degrees fibres]
k_overlap = 0.75          #(fig 5.5)
BOR = 0.016

rho_hydr = 70 # [kg/m3]

d_H_vap = 446100 #J/kg
BOR_percentage = 0.016 # [-]
d_T = 293.15 # temp diff between 20 K and 40 C
h_out = 30 #[W/m2K]

def multi_cell_v(m, n, p, R):
    """"Function to calculate multi cell volume in L."""
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


def multi_cell_m_incl_insulation(m,n,p, R, t_ins, Mass):
    """"Function to calculate multi cell mass in kg."""
    S, M_structural, t = multi_cell_s(m,n,p, R)

    volume = multi_cell_v(m, n, p, R )
    tanks_needed = math.ceil( Mass/(volume * rho_hydr))

    width = m * R * 2 - (m-1) * (y * R)
    length = n * R * 2 - (n-1) * (y * R)
    height = p * R * 2 - (p-1) * (y * R)

    dimensions = [width,length,height]

    thick = True
    t_ins = 0.01
    time = 0
    while thick:
        S_cryo, M_throwaway, t_throwaway = multi_cell_s(m,n,p, R - t_ins)
        Q = (Mass * (1+BOR_percentage) - Mass)  * d_H_vap
        U = Q / S_cryo / d_T
        t_ins = (1 / U - t / k_comp + 1 / h_out) * k_ins
        time = time +1

        if time > 10:
            thick = False
            print("t_ins",t_ins)

    M_insulation = rho_ins * t_ins * S_cryo
    M_total = M_structural + M_insulation
    return M_total, M_structural, M_insulation, S, S_cryo, t_ins, tanks_needed, dimensions




class Material:
    def __init__(self, K, E, rho):
        self.K = K      # [MPa]
        self.E = E      # [GPa]
        self.rho = rho  # [kg/m^3]

aluminium_2219 = Material(172.4*10**6, 73.8 * 10 ** 9, 2825)

class Tank:
    def __init__(self, m, n, p, R, t_ins, Mass):
        self.mass_total = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[0]
        self.mass_tank = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[1]
        self.mass_insulation = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[2]
        self.surface_outside = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[3]
        self.surface_inside_insulation = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[4]
        self.volume = multi_cell_v(m, n, p, R - multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[5]) * 1000
        self.number = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[6]
        self.dimensions = multi_cell_m_incl_insulation(m, n, p, R, t_ins, Mass)[7]



a  = Tank(8, 4, 2, 0.2, 0 , 100)

print('v:', a.volume)
print('m_tank:', a.mass_tank)
print('m_ins:', a.mass_insulation)
print('m_tot:', a.mass_total)
print('number of tanks:', a.number)
print('v_total_all:', a.number*a.volume)
print('M_total_all:', a.number*a.mass_total)
print('Dimensions:', a.dimensions)


#b = multi_cell_v(2, 2, 1, 0.145, 0.200, 0.029, 0.06047, 0.75, 0.00397)

# 0.037







