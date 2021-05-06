from math import *
import numpy as np

#Optimum lift coefficient, Cl_opt
def calculate_cl_opt(cd_0, A, e):
    """Calculates optimum lift coefficient, Cl_opt. Inputs are zero lift drag coefficient Cd_0, aspect ratio A
    and Oswald efficiency factor e."""
    return sqrt((1 / 3) * cd_0 * A * e)

#Cruise speed
def calculate_v_cruise(W, S, rho, cl_opt):
    """Calculates the cruise velocity, V_cruise. Inputs are weight W, wing surface area S, air density rho
    and optimum lift coefficient Cl_opt."""
    return sqrt((2 * W / S * rho * cl_opt))

#Required power
def calculate_p_cruise(W, cd_0, rho, v_cruise, S, A, e):
    """Calculates the required power, P_cruise, to fly at cruise speed. Inputs are weight W, zero lift drag
    coefficient Cd_0, air density rho, the cruise velocity v_cruise, wing surface area S, aspect ratio A and
    Oswald efficiency factor e."""
    return W * (((cd_0 * (1 / 2) * rho * v_cruise**3) / (W / S)) + ((W / S) * (1 / pi * A * e * (1 / 2) * rho *
                                                                               v_cruise)))
#Range
def calculate_range(specific_energy, m_energy, m, L, D, efficiency_total):
    """Calculates the range, R. Inputs are specific energy of the energy source E^*, mass of the energy source m_energy,
    mass of the aircraft m, lift L, drag D and total efficiency eta_total."""
    return specific_energy * (m_energy / m) * (1 / 9.81) * (L / D) * efficiency_total

#Climb rate
def calculate_max_climb_rate(p_max, W, S, cd_0, rho, A, e):
    """Calculates the maximum climb rate, c. Inputs are maximum power p_max, weight W, wing surface area S, zero lift
    drag coefficient Cd_0, air density rho, aspect ratio A and Oswald efficiency factor e."""
    cl = sqrt(3 * cd_0 * pi * A * e)
    cd = 4 * cd_0
    return (p_max / W) - ((sqrt(W / S) * sqrt(2)) / ((cl**(3 / 2) / cd) * sqrt(rho)))

#Climb gradient
def calculate_climb_gradient(c, v):
    """Calculates the climb gradient, G. Inputs are climb rate c and velocity v."""
    return c / v

#Runway length take-off
def calculate_runway_length_takeoff(v_final, v_initial, T, W, mu, rho, S, cl_takeoff, cd_0, A, e):
    """Calculates the runway length for take-off, s_g. Inputs are final velocity v_final, initial velocity v_initial,
    thrust T, weight W, ground friction constant mu, air density rho, wing surface area S, lift coefficient at take-off
    cl_takeoff, zero lift drag coefficient cd_0, aspect ratio A and Oswald efficiency factor e."""
    k = 1 / (pi * A * e)
    k_t = (T / W) - mu
    k_a = (rho / (2 * (W / S))) * (mu * cl_takeoff - cd_0 - k * cl_takeoff**2)
    safety_margin = 1.15
    return (1 / (2 * 9.81 * k_a)) * np.log(((k_t + k_a * v_final**2) / (k_t + k_a * v_initial**2))) * safety_margin

#Take-off parameter
def calculate_takeoff_parameter(W, S, p_max, cl_takeoff):
    """Calculates the take-off parameter, TOP_prop. Inputs are weight W, wing surface area S, maximum power p_max and
    lift coefficient at take-off cl_takeoff."""
    sigma = 1
    return (W / S) * (W / p_max) * (1 / cl_takeoff) * (1 / sigma)

#Runway length landing
def calculate_runway_length_landing(W, rho, A, cl_max):
    """Calculates the runway length for landing, s_l. Inputs are weight W, air density rho, aspect ratio A and maximum
    lift coefficient cl_max."""
    v_stall = sqrt((2 * W) / rho * A * cl_max)
    return 0.5915 * v_stall**2


