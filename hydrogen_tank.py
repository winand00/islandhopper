import numpy as np
eta_fc = 0.5
eta_total = 0.85 * eta_fc * 0.9 * 0.9

LD = 8.97
E_h2 = 120*10**6
payload = 1814.06
mass_max_payload = 8618
mass_no_payload = 8618 - payload
max_range = 1.852 * 300 *10**3
payload_range = 1.852 * 200 *10**3
rho_liquid = 70
rho_pressure = 40
energy_margin = 1.3
eta_storage_liquid = 0.2
eta_storage_pressure = 0.1
max_power = 1383.229 * 10**3
cruise_power = 1244.048 * 10**3
rho_pmad = 10000
rho_em = 5000
rho_fc = 3400
rho_comp = 2000
rho_cool = 2000
T_air = 298.15
T_h2 = 20
dT = T_air - T_h2


def diameter(volume):
    return 2*np.cbrt(volume/(4/3*np.pi))


def cylinder_length(volume):
    D_fuselage = 1.8
    return volume/(np.pi*(D_fuselage/2)**2)


def cooling_power():
    P_heat_rejected = (1/eta_fc-1) * max_power
    f_dT = 0.0038 * (T_air/dT)**2 + 0.0352 * T_air/dT + 0.1817
    return (0.371 * P_heat_rejected + 1.33) * f_dT


def tank_sizing(range, mass, rho_h2, eta_storage, rho_type, cooling):
    g = 9.81
    cooling_fraction = 1 + cooling_power()/max_power
    maximum_power = max_power
    mass_h2 = energy_margin * range/(eta_total*LD*E_h2)*mass*g
    if cooling:
        mass_h2 *= cooling_fraction
        maximum_power += cooling_power()
    volume_tank = mass_h2 / (rho_h2 * 0.5)
    length_tank = cylinder_length(volume_tank)
    mass_tank_hydrogen = mass_h2 / eta_storage
    mass_other_systems = maximum_power/rho_fc + maximum_power/rho_pmad + maximum_power/rho_em + cooling_power()/rho_type
    mass_total = mass_other_systems + mass_tank_hydrogen
    return mass_h2, mass_tank_hydrogen, volume_tank, length_tank, mass_total


print("Pressure, max payload", tank_sizing(payload_range, mass_max_payload, rho_pressure, eta_storage_pressure, rho_comp, False))
print("Liquid, max payload", tank_sizing(payload_range, mass_max_payload, rho_liquid, eta_storage_liquid, rho_cool, True))

print("Pressure, max range", tank_sizing(max_range, mass_no_payload, rho_pressure, eta_storage_pressure, rho_comp, False))
print("Liquid, max range", tank_sizing(max_range, mass_no_payload, rho_liquid, eta_storage_liquid, rho_cool, True))