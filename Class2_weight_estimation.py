from math import *
m_to_ft = 3.2808399
m2_to_ft2 = 10.7639104
kg_to_lbs = 2.20462262



def get_weight_wing(W_dg,N_z,S_w,A,t_c_root,Lambda,Sweep_angle,S_csw):
    S_w = S_w * m2_to_ft2
    W_dg = W_dg * kg_to_lbs
    S_csw = S_csw * m2_to_ft2
    W_wing = 0.0051 * (W_dg*N_z)**(0.557) * S_w**(0.649) * A**(0.5) * (t_c_root)**(-0.4) * (1 + Lambda)**(0.1) * cos(Sweep_angle)**(-1) * S_csw**0.1
    return W_wing/kg_to_lbs

def get_weight_avionics(W_uav):
    W_uav = W_uav * kg_to_lbs
    return (1.73 * W_uav**(0.983))/kg_to_lbs

def get_weight_landing_gear(K_mp,N_l,W_l,L_m,L_n,N_mw,N_nw,N_mss,V_stall):
    W_l = W_l * kg_to_lbs
    L_m = L_m * m_to_ft
    L_n = L_n * m_to_ft
    W_mlg = 0.095*(N_l*W_l)**(0.768) * (L_m/12)**(0.409)
    W_nlg = 0.125 * (N_l * W_l) ** (0.566) * (L_n/12) ** (0.845)
    return (W_mlg + W_nlg)/kg_to_lbs

#def get_weight_tail():




# Inputs
#Wing inputs:
S_w = 45             #Wing surface [m^2]
W_dg = 8618             #Design gross weight [kg]        (took 12% of MTOW for now)
N_z = 4.3913             #Ultimate load factor (1.5* limit load factor)
A = 9                #Aspect ratio
t_c_root =   0.16       #Thickness to chord ratio at root [m]
Lambda =  0.5          #Wing taper ratio
Sweep_angle = 0      #Sweep angle at 25% MAC
S_csw = 2.54           #Control surface area [m^2]

#Avionics inputs
W_uav = 1   #Uninstalled avionics weight

#Landing gear inputs
K_mp
L_n =    1       #nose gear length [m]
L_m =    1       #Length of main landing gear
N_l =  4.39         #Ultimate landing load landing factor       (Used N_Z for now)
W_l = 0.05*8618          #Landing design gross weight [kg]    (took 5% of MTOW for now)

weight_wing = get_weight_wing(W_dg,N_z,S_w,A,t_c_root,Lambda,Sweep_angle,S_csw)
weight_avionics = get_weight_avionics(W_uav)
weight_landing_gear = get_weight_landing_gear(N_l,W_l,L_m,L_n)
print("Weight wing =", weight_wing)
print("Weight avionics =", weight_avionics)
print("Weight Landing gear=", weight_landing_gear)