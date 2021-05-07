from math import *
import numpy as np


# Optimum lift coefficient, Cl_opt
def calculate_cl_opt(cd_0, A, e):
    """Calculates optimum lift coefficient, Cl_opt. Inputs are zero lift drag coefficient Cd_0, aspect ratio A
    and Oswald efficiency factor e."""
    return sqrt((1 / 3) * cd_0 * A * e)


# Cruise speed
def calculate_v_cruise(W, S, rho, cl_opt):
    """Calculates the cruise velocity, V_cruise. Inputs are weight W, wing surface area S, air density rho
    and optimum lift coefficient Cl_opt."""
    return sqrt((2 * W / (S * rho * cl_opt)))

# Required power
def calculate_p_cruise(W, cd_0, rho, v_cruise, S, A, e):
    """Calculates the required power, P_cruise, to fly at cruise speed. Inputs are weight W, zero lift drag
    coefficient Cd_0, air density rho, the cruise velocity v_cruise, wing surface area S, aspect ratio A and
    Oswald efficiency factor e."""
    return W * (((cd_0 * (1 / 2) * rho * v_cruise ** 3) / (W / S)) + ((W / S) * (1 / (pi * A * e * (1 / 2) * rho *
                                                                                 v_cruise))))

# Range
def calculate_range(specific_energy, m_energy, m, L_over_D, efficiency_total):
    """Calculates the range, R. Inputs are specific energy of the energy source E^*, mass of the energy source m_energy,
    mass of the aircraft m, lift L, drag D and total efficiency eta_total."""
    return specific_energy * (m_energy / m) * (1 / 9.81) * L_over_D * efficiency_total


# Climb rate
def calculate_max_climb_rate_and_gradient(p_max, W, S, cd_0, rho, A, e):
    """Calculates the maximum climb rate c, the climb rate G at maximum climb rate and the velocity V at maximum climb
    rate. Inputs are maximum power p_max, weight W, wing surface area S, zero lift drag coefficient Cd_0, air density
    rho, aspect ratio A and Oswald efficiency factor e."""
    cl = sqrt(3 * cd_0 * pi * A * e)
    cd = 4 * cd_0
    max_climb_rate = (p_max / W) - ((sqrt(W / S) * sqrt(2)) / ((cl ** (3 / 2) / cd) * sqrt(rho)))
    max_climb_gradient = (p_max / W) * (1 / sqrt((W / S) * (2 / rho) * (1 / cl))) - cd / cl
    velocity_at_max_climb_rate = max_climb_rate / max_climb_gradient
    return max_climb_rate, max_climb_gradient, velocity_at_max_climb_rate


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


# Take-off parameter
def calculate_takeoff_parameter(W, S, p_max, cl_takeoff):
    """Calculates the take-off parameter, TOP_prop. Inputs are weight W, wing surface area S, maximum power p_max and
    lift coefficient at take-off cl_takeoff."""
    sigma = 1
    return (W / S) * (W / p_max) * (1 / cl_takeoff) * (1 / sigma)


# Take-off distance
def calculate_takeoff_distance(TOP):
    return 0.0577 * TOP ** 2 + 8.6726 * TOP


# Runway length landing
def calculate_landing_distance(W, rho_sealevel, S, cl_max):
    """Calculates the runway length for landing, s_l. Inputs are weight W, air density rho, aspect ratio A and maximum
    lift coefficient cl_max."""
    return W * 2 * 0.5915 / (cl_max * rho_sealevel * S)

# Refuel time
def calculate_refuel_time(tank_capacity, fill_speed):
    """Calculates refuel time t_rf. Inputs are tank capacity and fill speed."""
    return tank_capacity / fill_speed


def calculate_charging_time(E_batt, p_ch, efficiency_ch):
    """Calculates charging time t_ch. Inputs are battery capacity E_batt, charging power p_ch and charging efficiency
    eta_ch."""
    return E_batt / (p_ch * efficiency_ch)


# Noise
def calculate_max_sound_pressure_level(p_br, D, M_t, B, N, r):
    """Calculates the maximum sound pressure level, SPL_max. Inputs are break shaft power p_br, propeller diameter D,
    rotational tip Mach number M_t, number of blades per propeller B, number of propellers N and distance to the
    propeller r."""
    return 83.4 + 15.3 * np.log10(p_br) - 20 * np.log10(D) + 38.5 * M_t - 3 * (B - 2) + 10 * np.log10(N) - 20 * \
           np.log10(r)


# Design efficiency
def calculate_design_efficiency(E_d, efficiency_pt, efficiency_r, E_T):
    """Calculates the design efficiency GI_1. Inputs are the amount of energy per kg energy source material E_d, total
    efficiency of the powertrain eta_pt, percentage of the total used energy that is recovered for other systems eta_r
    and the total amount of thrust energy T."""
    return (E_d * efficiency_pt * (1 + efficiency_r)) / E_T


# Tool
def tool(cd_0, A, e, W, rho, S, specific_energy, m_energy, m, L_over_D, efficiency_total, p_max, cl_takeoff, cl_max,
         p_br, D, M_t, B, N, r, E_d, efficiency_pt, efficiency_r, E_T, battery=True):
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
    p_cruise = calculate_p_cruise(W, cd_0, rho, v_cruise, S, A, e)
    max_range = calculate_range(specific_energy, m_energy, m, L_over_D, efficiency_total)

    print(f"**************CRUISE CHARACTERISTICS******************** \n"
          f"Optimum lift coefficient Cl_opt = {round(cl_opt, 2)} [-]\n"
          f"Cruise speed                    = {round(v_cruise, 2)} [m/s] \n"
          f"Required cruise power           = {round(p_cruise, 2)} [Watt] \n"
          f"Max range                       = {round(max_range, 2)} [m]")

    max_climb_rate, max_climb_gradient, climb_velocity = calculate_max_climb_rate_and_gradient(p_max, W, S, cd_0, rho,
                                                                                               A, e)
    #runway_length_takeoff = calculate_runway_length_takeoff(v_final, v_initial, T, W, mu, rho, S, cl_takeoff, cd_0, A, e)
    TOP = calculate_takeoff_parameter(W, S, p_max, cl_takeoff)
    runway_length_landing = calculate_runway_length_landing(W, rho, A, cl_max)

    print(f"***********TAKE-OFF/LANDING CHARACTERISTICS************* \n"
          f"Max climb rate                  = {round(max_climb_rate, 2)} [m/s] \n"
          f"Climb gradient @max climb rate  = {round(degrees(atan(max_climb_gradient)), 2)} [deg] \n"
          f"Velocity       @max climb rate  = {round(climb_velocity, 2)} [m/s] \n"
          f"Take-Off Parameter (TOP)        = {round(TOP, 2)} [N^2/Wm^2] \n"
          #f"Runway length @take-off         = {round(runway_length_takeoff, 2)} [m] \n"
          f"Runway length @landing          = {round(runway_length_landing, 2)} [m]")

    if battery:
        E_batt = float(input('E_batt [J]    = '))
        p_ch = float(input('p_ch [W]     = '))
        efficiency_ch = float(input('efficiency_ch = '))
        charge_time = calculate_charging_time(E_batt, p_ch, efficiency_ch)
        print(f"***********OTHER CHARACTERISTICS************* \n"
              f"Battery charge time             = {round(charge_time, 2)} [sec]")

    else:
        tank_capacity = float(input('tank capacity [kg] = '))
        fill_speed = float(input('fill speed  [kg/s] = '))
        refuel_time = calculate_refuel_time(tank_capacity, fill_speed)
        print(f"***********OTHER CHARACTERISTICS************* \n"
              f"Refuel time                     = {round(refuel_time, 2)} [sec]")

    max_SPL = calculate_max_sound_pressure_level(p_br, D, M_t, B, N, r)
    design_efficiency = calculate_design_efficiency(E_d, efficiency_pt, efficiency_r, E_T)

    print(f"Maximum sound pressure level    = {round(max_SPL, 2)} [dB] \n"
          f"Design efficiency parameter GI1 = {round(design_efficiency, 2)} [-]")


# Inputs jetstream 31
if __name__ == "__main__":
    tool(0.0376, 10, 0.75, 68179.5, 0.54895, 25.08, 46200000, 1483, 6950, 8.16, 0.65, 1074000, 1.8, 2, 771, 2.69,
         0.9, 4, 2, 100, 46200000, 0.7, 0.3, 4567890, False)