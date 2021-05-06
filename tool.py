from math import *

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


