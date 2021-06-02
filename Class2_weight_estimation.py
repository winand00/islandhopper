from math import *
m_to_ft = 3.2808399
m2_to_ft2 = 10.7639104
kg_to_lbs = 2.20462262
m_to_in = 39.3700787



def get_weight_wing(W_dg,N_z,S_w,A,t_c_root,Lambda,Sweep_angle,S_csw):
    S_w = S_w * m2_to_ft2
    W_dg = W_dg * kg_to_lbs
    S_csw = S_csw * m2_to_ft2
    W_wing = 0.0051 * (W_dg*N_z)**(0.557) * S_w**(0.649) * A**(0.5) * (t_c_root)**(-0.4) * (1 + Lambda)**(0.1) * cos(Sweep_angle)**(-1) * S_csw**0.1
    return W_wing/kg_to_lbs

def get_weight_avionics(W_uav):
    W_uav = W_uav * kg_to_lbs
    return (1.73 * W_uav**(0.983))/kg_to_lbs

def get_weight_landing_gear(K_mp,K_np,N_l,W_l,L_m,L_n,N_mw,N_nw,N_mss,V_stall):
    W_l = W_l * kg_to_lbs
    L_m = L_m * m_to_in
    L_n = L_n * m_to_in
    V_stall = V_stall * m_to_ft
    W_mlg = 0.0106 * K_mp*W_l**(0.888) * N_l**(0.25) * L_m**(0.4) * N_mw**(0.321) * N_mss**(-0.5) * V_stall**(0.1)
    W_nlg = 0.032*K_np*W_l**(0.646) *N_l**(0.2) * L_n**(0.5)*N_nw**(0.45)
    return (W_mlg + W_nlg)/kg_to_lbs

def get_weight_fuselage(K_door,K_lg,W_dg,N_z,L,S_f,K_ws,L_over_D):
    K_ws = 0.75 * ((1+2*Lambda)/(1+Lambda)) * B_w * tan(Sweep_angle/L_fus)

#Fuselage inputs
K_door = 1     #1.0 if no cargo door; = 1.06 if one side cargo door; = 1.12 if two side cargo doors; = 1.12 if aft clamshell door; = 1.25 if two side cargo doors and aft clamshell door
K_lg = 1.12    # 1.12 if fuselage-mounted main landing gear;= 1.0 otherwise
L_fus =         #Fuselage length
S_f =           #Fuselage wetted area
L_over_D = 13   #

# Inputs
W_dg = 8618             #Design gross weight [kg]
N_z = 4.3913             #Ultimate load factor (1.5* limit load factor)

#Wing inputs:
S_w = 45             #Wing surface [m^2]
A = 9                #Aspect ratio
t_c_root =   0.16       #Thickness to chord ratio at root [m]
Lambda =  0.5          #Wing taper ratio
Sweep_angle = 0      #Sweep angle at 25% MAC
S_csw = 2.54           #Control surface area [m^2]

#Avionics inputs
W_uav = 1   #Uninstalled avionics weight

#Landing gear inputs
K_mp = 1         #1.126 for kneeling gear;= 1.0 otherwise
K_np = 1         #1.15 for kneeling gear;= 1.0 otherwise
N_mw =  2        #Number of main wheels
N_nw =  1        #number of nose wheels
N_mss = 2        #Number of main gear shock struts   (No clue how much)
V_stall = 58     #Stall velocity
L_n =   2       #nose gear length [m]
L_m =    2       #Length of main landing gear [m]
N_l = 4.3913         #Ultimate landing load landing factor       (Used N_Z for now)
W_l = 8618          #Landing design gross weight [kg]    (took 5% of MTOW for now)

weight_wing = get_weight_wing(W_dg,N_z,S_w,A,t_c_root,Lambda,Sweep_angle,S_csw)
weight_avionics = get_weight_avionics(W_uav)
weight_landing_gear = get_weight_landing_gear(K_mp,K_np,N_l,W_l,L_m,L_n,N_mw,N_nw,N_mss,V_stall)
print("Weight wing =", weight_wing)
print("Weight avionics =", weight_avionics)
print("Weight Landing gear=", weight_landing_gear)