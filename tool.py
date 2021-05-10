from math import *
import numpy as np

"""
On the bottom of the file you can put in your own values for all parameters
"""


# Optimum lift coefficient, Cl_opt#
def calculate_cl_opt(cd_0, A, e):
    """Calculates optimum lift coefficient, Cl_opt. Inputs are zero lift drag coefficient Cd_0, aspect ratio A
    and Oswald efficiency factor e."""
    return sqrt((1 / 3) * cd_0 * A * e)


# Cruise speed#
def calculate_v_cruise(W, S, rho, cl_opt):
    """Calculates the cruise velocity, V_cruise. Inputs are weight W, wing surface area S, air density rho
    and optimum lift coefficient Cl_opt."""
    return sqrt((2 * W / (S * rho * cl_opt)))

# Required power
#def calculate_p_cruise(W, cd_0, rho, v_cruise, S, A, e):
#    """Calculates the required power, P_cruise, to fly at cruise speed. Inputs are weight W, zero lift drag
#    coefficient Cd_0, air density rho, the cruise velocity v_cruise, wing surface area S, aspect ratio A and
#    Oswald efficiency factor e."""
#    return W * (((cd_0 * (1 / 2) * rho * v_cruise ** 3) / (W / S)) + ((W / S) * (1 / (pi * A * e * (1 / 2) * rho *
#                                                                                 v_cruise))))

# Range#
def calculate_range(specific_energy, m_energy, m, L_over_D, efficiency_fuelcell, efficiency_prop):
    """Calculates the range, R. Inputs are specific energy of the energy source E^*, mass of the energy source m_energy,
    mass of the aircraft m, lift L, drag D and total efficiency eta_total."""
    return specific_energy * (m_energy / m) * (1 / 9.81) * L_over_D * efficiency_fuelcell * efficiency_prop


# Climb rate#
def calculate_max_climb_rate_and_gradient(p_max, W, S, cd_0, rho, A, e,efficiency_prop):
    """Calculates the maximum climb rate c, the climb rate G at maximum climb rate and the velocity V at maximum climb
    rate. Inputs are maximum power p_max, weight W, wing surface area S, zero lift drag coefficient Cd_0, air density
    rho, aspect ratio A and Oswald efficiency factor e."""
    cl = sqrt(3 * cd_0 * pi * A * e)
    cd = 4 * cd_0
    max_climb_rate = efficiency_prop*(p_max / W) - ((sqrt(W / S) * sqrt(2)) / ((cl ** (3 / 2) / cd) * sqrt(rho)))
    max_climb_gradient = efficiency_prop*(p_max / W) * (1 / sqrt((W / S) * (2 / rho) * (1 / cl))) - (cd / cl)
    #velocity_at_max_climb_rate = max_climb_rate / max_climb_gradient
    return max_climb_rate, max_climb_gradient


# Runway length take-off
#def calculate_runway_length_takeoff(v_final, v_initial, T, W, mu, rho, S, cl_takeoff, cd_0, A, e):
#    """Calculates the runway length for take-off, s_g. Inputs are final velocity v_final, initial velocity v_initial,
#    thrust T, weight W, ground friction constant mu, air density rho, wing surface area S, lift coefficient at take-off
#    cl_takeoff, zero lift drag coefficient cd_0, aspect ratio A and Oswald efficiency factor e."""
#    k = 1 / (pi * A * e)
#    k_t = (T / W) - mu
#    k_a = (rho / (2 * (W / S))) * (mu * cl_takeoff - cd_0 - k * cl_takeoff ** 2)
#    safety_margin = 1.15
#    return (1 / (2 * 9.81 * k_a)) * np.log(((k_t + k_a * v_final ** 2) / (k_t + k_a * v_initial ** 2))) * safety_margin


# Take-off parameter#
def calculate_takeoff_parameter(W, S, p_max, cl_takeoff):
    """Calculates the take-off parameter, TOP_prop. Inputs are weight W, wing surface area S, maximum power p_max and
    lift coefficient at take-off cl_takeoff."""
    sigma = 1
    return (W / S) * (W / p_max) * (1 / cl_takeoff) * (1 / sigma)


# Take-off distance#
def calculate_takeoff_distance(TOP):
    return 0.0577 * TOP ** 2 + 8.6726 * TOP


# Runway length landing#
def calculate_landing_distance(W, rho_sealevel, S, cl_max):
    """Calculates the runway length for landing, s_l. Inputs are weight W, air density rho, aspect ratio A and maximum
    lift coefficient cl_max."""
    return W * 2 * 0.5915 / (cl_max * rho_sealevel * S)

# Refuel time#
def calculate_refuel_time(tank_capacity, fill_speed):
    """Calculates refuel time t_rf. Inputs are tank capacity and fill speed."""
    return tank_capacity / fill_speed


def calculate_charging_time(E_batt, p_ch, efficiency_ch):
    """Calculates charging time t_ch. Inputs are battery capacity E_batt, charging power p_ch and charging efficiency
    eta_ch."""
    return E_batt / (p_ch * efficiency_ch)


# Noise#
def calculate_max_sound_pressure_level(D, B, N, p_max, efficiency_prop):
    """Calculates the maximum sound pressure level, SPL_max. Inputs are break shaft power p_br, propeller diameter D,
    rotational tip Mach number M_t, number of blades per propeller B, number of propellers N and distance to the
    propeller r."""
    return 83.4 + 15.3 * np.log10(p_max / (N * efficiency_prop)) - 20 * np.log10(D) + 38.5 * 0.9 - 3 * (B - 2) + 10 * np.log10(N) - 20 * \
           np.log10(100)


# Design efficiency#
def calculate_design_efficiency(specific_energy, efficiency_fuelcell, efficiency_prop, efficiency_r, m_energy):
    """Calculates the design efficiency GI_1. Inputs are the amount of energy per kg energy source material E_d, total
    efficiency of the powertrain eta_pt, percentage of the total used energy that is recovered for other systems eta_r
    and the total amount of thrust energy T."""
    return (specific_energy * efficiency_fuelcell * efficiency_prop * (1 + efficiency_r)) / (specific_energy * m_energy)


# Tool
def tool(cd_0, A, e, W, rho, rho_sealevel, S, specific_energy, m_energy, m, L_over_D, efficiency_fuelcell,
         efficiency_prop, p_max, cl_takeoff, cl_max, D, B, N, efficiency_r, battery=True):
    # (cd_0, A, e, W, rho, S, specific_energy, m_energy, m, L_over_D, efficiency_total, p_max, v_initial, v_final, T,
    #  mu, cl_takeoff, cl_max, p_br, D, M_t, B, N, r, E_d, efficiency_pt, efficiency_r, E_T, battery=True)

    # print("Enter parameters")
    # cd_0 = float(input("Zero lift drag, Cd_0:"))
    # W = float(input("Weight, W:"))
    # A = float(input("Aspect ratio, A:"))
    # e = float(input("Oswald efficiency factor, e:"))
    # rho = float(input("Air density, rho:"))
    # S = float(input("Wing surface area, S:"))
    # specific_energy = float(input("Specific energy, E^*:"))
    # m_energy = float(input("Mass energy, m_energy:"))
    # m = float(input("Mass, m:"))
    # L_over_D = float(input("Lift over drag, L/D:"))
    # efficiency_total = float(input("Total efficiency, eff_tot:"))
    # p_max = float(input("Maximum power (includes efficiency skim), P_max:"))
    # v_initial = float(input("Take-off initial velocity, v_i:"))
    # v_final = float(input("Take-off final velocity, v_f:"))
    # T = float(input("Thrust, T:"))
    # mu = float(input("Ground friction constant, mu:"))
    # cl_takeoff = float(input("Lift coefficient @take-off, Cl_takeoff:"))
    # cl_max = float(input("Maximum lift coefficient, Cl_max"))
    # p_br = float(input("Break shaft power, P_br:"))
    # D = float(input("Propeller diameter, D:"))
    # M_t = float(input("Rotational tip Mach number, M_t:"))
    # B = float(input("Number of blades per propeller, B:"))
    # N = float(input("Number of propellers, N:"))
    # r = float(input("Distance to the propeller, r:"))
    # E_d = float(input("Amount of energy per kg energy source material, E_d:"))
    # efficiency_pt = float(input("Efficiency of the powertrain, eta_pt:"))
    # efficiency_r = float(input("Percentage of the total used energy that is recovered for other systems, eta_r:"))
    # E_T = float(input("Total amount of thrust energy, E_T:"))
    # battery = input("Battery aircraft (True/False):")

    cl_opt = calculate_cl_opt(cd_0, A, e)
    v_cruise = calculate_v_cruise(W, S, rho, cl_opt)
    #p_cruise = calculate_p_cruise(W, cd_0, rho, v_cruise, S, A, e)
    max_range = calculate_range(specific_energy, m_energy, m, L_over_D, efficiency_fuelcell, efficiency_prop)

    print(f"**************CRUISE CHARACTERISTICS******************** \n"
          #f"Optimum lift coefficient Cl_opt = {round(cl_opt, 2)} [-]\n"
          f"Cruise speed                    = {round(v_cruise, 2)} [m/s] \n"
          #f"Required cruise power           = {round(p_cruise, 2)} [Watt] \n"
          f"Max range                       = {round(max_range, 2)} [m]")

    max_climb_rate, max_climb_gradient = calculate_max_climb_rate_and_gradient(p_max, W, S, cd_0, rho, A, e,efficiency_prop)
    #runway_length_takeoff = calculate_runway_length_takeoff(v_final, v_initial, T, W, mu, rho, S, cl_takeoff, cd_0, A, e)
    #TOP = calculate_takeoff_parameter(W, S, p_max, cl_takeoff)
    runway_length_landing = calculate_landing_distance(W, rho_sealevel, S, cl_max)
    TOP = calculate_takeoff_parameter(W, S, p_max, cl_takeoff)
    takeoff_distance = calculate_takeoff_distance(TOP)

    print(f"***********TAKE-OFF/LANDING CHARACTERISTICS************* \n"
          f"Max climb rate                  = {round(max_climb_rate, 2)} [m/s] \n"
          f"Climb gradient @max climb rate  = {round(degrees(atan(max_climb_gradient)), 2)} [deg] \n"
          #f"Velocity       @max climb rate  = {round(climb_velocity, 2)} [m/s] \n"
          f"Take-Off Parameter (TOP)        = {round(TOP, 2)} [N^2/Wm^2] \n"
          f"Take-off distance         = {round(takeoff_distance, 2)} [m] \n"
          #f"Runway length @take-off         = {round(runway_length_takeoff, 2)} [m] \n"
          f"Landing distance          = {round(runway_length_landing, 2)} [m]")

    if battery:
        E_batt = float(input('E_batt [J]    = '))
        p_ch = float(input('p_ch [W]     = '))
        efficiency_ch = float(input('efficiency_ch [-]     = '))
        charge_time = calculate_charging_time(E_batt, p_ch, efficiency_ch)
        print(f"***********OTHER CHARACTERISTICS************* \n"
              f"Battery charge time             = {round(charge_time, 2)} [sec]")

    else:
        tank_capacity = float(input('tank capacity [kg] = '))
        fill_speed = float(input('fill speed  [kg/s] = '))
        refuel_time = calculate_refuel_time(tank_capacity, fill_speed)
        print(f"***********OTHER CHARACTERISTICS************* \n"
              f"Refuel time                     = {round(refuel_time, 2)} [sec]")

    max_SPL = calculate_max_sound_pressure_level(D, B, N, p_max, efficiency_prop)
    design_efficiency = calculate_design_efficiency(specific_energy, efficiency_fuelcell, efficiency_prop,
                                                    efficiency_r, m_energy)

    print(f"Maximum sound pressure level    = {round(max_SPL, 2)} [dB] \n"
          f"Design efficiency parameter GI1 = {round(design_efficiency, 2)} [-]")


# Inputs jetstream 31

cd_0 = 0.0376               # Zero lift drag
A = 10                      # Aspect ratio
e = 0.75                    # Oswald efficiency
W = 68179.5                 # Weight of aircraft [N]
rho = 0.54895               # Density at cruise h [kg/m3]
rho_sealevel = 1.225        # Density at sealevel [kg/m3]
S = 25.08                   # Wing  [m2]
specific_energy = 46200000  # Specific energy of fuel [J/kg]
m_energy = 1483             # Mass of fuel [kg]
m = 6950                    # Mass of aircraft [kg]
L_over_D = 8.16             # Lift over drag ratio
efficiency_fuelcell = 0.9   # Efficiency fuel cell
efficiency_prop = 0.8       # Efficiency propeller
p_max = 1074000             # Max power [W]
cl_takeoff = 1.8            # Lift coefficient at take-off
cl_max = 2                  # Maximum lift coefficient
D = 2.69                    # Propeller diameter
B = 4                       # Number of blades per propeller
N = 2                       # Number of propellers
efficiency_r = 0.1          # Fraction of total used energy that is recovered for other systems
battery = False             # Aircraft on batteries

if __name__ == "__main__":
   tool(cd_0, A, e, W, rho, rho_sealevel, S, specific_energy, m_energy, m, L_over_D, efficiency_fuelcell,
        efficiency_prop, p_max, cl_takeoff, cl_max, D, B, N, efficiency_r, battery)
