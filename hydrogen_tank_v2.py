import numpy as np

eta_fc = 0.5
eta_elecsys = 0.9
eta_em = 0.9
eta_pmad = 0.9
eta_total = 0.85 * eta_fc * eta_em * eta_pmad * eta_elecsys

LD = 12.041452736768791
E_h2 = 120 * 10 ** 6
payload = 1814.06
mass_max_payload = 8618
mass_no_payload = 8618 - payload
max_range = 1.852 * 300 * 10 ** 3 + 50.32 * 45 * 60
payload_range = 1.852 * 200 * 10 ** 3 + 50.32 * 45 * 60
rho_liquid = 70
rho_pressure = 40
energy_margin = 1
eta_storage_liquid = 0.2
eta_storage_pressure = 0.1
max_power = 1294262.6121018375
cruise_power = 1114089.2877240134
rho_pmad = 10000
rho_em = 5000
rho_fc = 3400
rho_comp = 2000
rho_cool = 2000
T_air = 298.15
T_h2 = 20
dT = T_air - T_h2

proof_pressure = 4 * 10**5
safety_factor = 1.5

a = c = 0.9


class Material:
    def __init__(self, K, E, rho):
        self.K = K      # [MPa]
        self.E = E      # [GPa]
        self.rho = rho  # [kg/m^3]


aluminium_2219 = Material(172.4*10**6, 73.8 * 10 ** 9, 2825)

# 6061-T6
#
# 2024-T4
# 2219-T87
#
# 5052-H38
# 5083-H38

#https://www.gasparini.com/en/blog/metals-and-materials-for-low-temperatures/#:~:text=These%20are%20mainly%20quenched%20and,be%20used%20at%20these%20temperatures.

'''
determine outer diameter fuselage:
include: tank thickness + insulation

calulate length:
with hydrogen mass and predetermined shape

Calculate wall thickness

Calculate tank mass




'''


def diameter(volume):
    return 2 * np.cbrt(volume / (4 / 3 * np.pi))


def cylinder_length(volume):
    D_fuselage = 1.8
    return volume / (np.pi * (D_fuselage / 2) ** 2)


def cooling_power():
    P_heat_rejected = (1 / eta_fc - 1) * max_power
    f_dT = 0.0038 * (T_air / dT) ** 2 + 0.0352 * T_air / dT + 0.1817
    return (0.371 * P_heat_rejected + 1.33) * f_dT


def tank_sizing(range, mass, rho_h2, eta_storage, rho_type, cooling):
    g = 9.81
    cooling_fraction = 1 + cooling_power() / max_power
    maximum_power = max_power
    mass_h2 = energy_margin * range / (eta_total * LD * E_h2) * mass * g
    mass_h2 *= cooling_fraction
    maximum_power += cooling_power
    return(mass_h2, maximum_power)

    '''
    volume_tank = mass_h2 / (rho_h2 * 0.5)
    length_tank = cylinder_length(volume_tank)
    mass_tank_hydrogen = mass_h2 / eta_storage
    mass_other_systems = maximum_power / rho_fc + maximum_power / rho_pmad + maximum_power / rho_em + cooling_power() / rho_type
    mass_total = mass_other_systems + mass_tank_hydrogen
    return mass_h2, mass_tank_hydrogen, volume_tank, length_tank, mass_total
    '''

def thickness(material):


    rightside = proof_pressure * ((a + c) / (2 * s_w)) * (1 + 2 * (1 + 3.6 * (proof_pressure /
                                              material.E) * ((a + c) / 2 * s_w)**3)*((a - c) / (a + c))+ (1 / 2))

    leftside = material.K/safety_factor

    pass




def multi_cell(m,n,p, R, d, R_fillet, s, k_overlap, t_junction):
    # m,n,p are number of cells in 3 directions
    # R is desired radius
    # d = y*R, choose y value
    # R_fillet = choose value, !!r = can also be used if found
    # !!s, cant find it
    # k_overlap = 0.75 #Used in thesis
    # t_junction is dependent on y

    P = 5.7 #[bar] ? figure 5-11
    # sigma_allowable = ?  They use reference [4] and [5] to find it
    t_ply = 0.00014 #[m] #Table 4-3

    dH_vap = 446.1 * 1000 # [J/kg] https://link.springer.com/referenceworkentry/10.1007%2F978-90-481-2642-2_327#:~:text=The%20heat%20of%20vaporization%20of,boiling%20point%20under%20standard%20pressure.
    h_out = 30 #[W/m2 K] p.117
    k_comp = 1
    k_ins = 1


    N_cells = m*n*p
    V_spheres = N_cells * 4/3 * np.pi * R**3

    N_junctions = (n*(m-1) + m*(n-1)) * p + m*n*(p-1)
    N_centers = ((m-1)*(n-1))*p + (p-1)*(n-1)

    V_lensjunction = (np.pi*(2*R-d)**2*(d**2+4*d*R))/(12*d)
    V_lenses = N_junctions*V_lensjunction

    #R_fillet = R*r

    theta_1 = np.arcsin((d/2)/(R+R_fillet)) # e1 4.42
    h_inter = R_fillet*np.sin(theta_1)
    R_ring = np.sqrt(R**2 - (d/2)**2) #eq 4.41


    R_ringfinal = R_ring + ((R*np.cos(theta_1)-R_ring)-(h_inter*np.tan(theta_1/2)))
    V_centerlens = (np.pi*(2*R_ringfinal - d)**2*(d**2+4*d*R_ringfinal))/(12*d)
    V_centers = N_centers * V_centerlens

    h_center = (2*R - d*2**0.5)/2
    r_cyl = max((4*s*k_overlap)/2/np.pi , 2*t_junction, h_center)
    l_cyl = np.sqrt(R_ringfinal ** 2 - (d / 2 - r_cyl) ** 2)
    V_cylinder = np.pi * r_cyl**2 * l_cyl
    V_cylinders = N_centers * V_cylinder

    A_triangle = ((2*R_fillet*np.sin(theta_1))*(R*np.cos(theta_1)-R_ring))/2
    A_cap = R_fillet**2/2*(2*theta_1-np.sin(2*theta_1))
    V_fillet = (A_triangle - A_cap)*(2*np.pi*R_ring)
    V_fillets = N_junctions * V_fillet

    V = V_spheres + V_centers + V_fillets - V_lenses - V_cylinders

    S_spheres = N_cells*4*np.pi*R**2
    h_cap = R - np.sqrt(R**2 - R_ring**2)
    S_lensjunction = 2*(2*h_cap)
    S_lenses = N_junctions * S_lensjunction

    gamma_cyl = np.arcsin(l_cyl/(2*R_ringfinal))
    S_torus = 2*np.pi * s * R_ring*((2*np.pi-gamma_cyl)/(2*np.pi))
    S_fillets = N_junctions*S_torus
    S_spherefillet = 2*(2*np.pi*R*h_inter)
    S_spheresfillets = N_junctions * S_spherefillet

    S_center = 2*(2*np.pi*R*h_center)
    S_centers = N_centers * S_center

    S_cylinder = 2*np.pi*(r_cyl*l_cyl+r_cyl**2)
    S_cylinders = N_centers* S_cylinder

    S = S_spheres + S_fillets + S_cylinders + S_centers - S_lenses - S_spheresfillets

    t = max(P*R/(2*sigma_allowable), 6*t_ply)

    M_structural = rho_comp * t *(S_spheres + S_centers - S_lenses - S_spheresfillets)
                    + rho_comp*t_junction*(S_fillets+S_cylinders)

    M_insulation = rho_ins * t_ins * S_cryo #!!!!!!!!!!!!!!!!!!!!!!!!
    M_cryo = M_structural + M_insulation

    print('Volume =', V*1000, '[L]')
    print('Surface=', S, '[m^2]')

    return(V,S)

multi_cell(2,2,1, 0.145, 0.200, 0.029, 0.06047, 0.75, 0.00397)







# print("Pressure, max payload", tank_sizing(payload_range, mass_max_payload, rho_pressure, eta_storage_pressure, rho_comp, False))
# print("Liquid, max payload", tank_sizing(payload_range, mass_max_payload, rho_liquid, eta_storage_liquid, rho_cool, True))

# print("Pressure, max range", tank_sizing(max_range, mass_no_payload, rho_pressure, eta_storage_pressure, rho_comp, False))
# print("Liquid, max range", tank_sizing(max_range, mass_no_payload, rho_liquid, eta_storage_liquid, rho_cool, True))