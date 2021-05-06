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



