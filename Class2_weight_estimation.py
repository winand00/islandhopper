from math import *
m_to_ft = 3.2808399
m2_to_ft2 = 10.7639104
m3_to_ft3 = 35.3146667
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

def get_weight_fuselage(K_door,K_lg,W_dg,N_z,L_fus,S_f,L_over_D,B_w,Lambda):
    W_dg = W_dg * kg_to_lbs
    L_fus = L_fus * m_to_ft
    S_f = S_f * m2_to_ft2
    B_w = B_w * m_to_ft
    K_ws = 0.75 * ((1+2*Lambda)/(1+Lambda)) * B_w * tan(Sweep_angle/L_fus)
    W_fus = 0.3280*K_door*K_lg*(W_dg*N_z)**(0.5) * L_fus**(0.25) * S_f**(0.302) * (1+K_ws)**(0.04)*(L_over_D)**(0.1)
    return W_fus/kg_to_lbs


def get_weight_vert_tail(H_t,H_v,W_dg,N_z,L_t,S_vt,Sweep_angle_vt,A_v,t_c_root):
    W_dg = W_dg * kg_to_lbs
    L_t = L_t * m_to_ft
    H_t = H_t * m_to_ft
    H_v = H_v * m_to_ft
    S_vt = S_vt * m2_to_ft2
    K_z = L_t
    W_vt = 0.0026 * (1+H_t/H_v)**(0.225) * W_dg**(0.556) * N_z ** (0.536) * L_t**(-0.5) * S_vt**(0.5) * K_z**(0.875) * (cos(Sweep_angle_vt))**(-1) * A_v**(0.35) * t_c_root**(-0.5)
    return W_vt/kg_to_lbs

def get_weight_hor_tail(K_uht,F_w,B_h,W_dg,N_z,S_ht,L_t,Sweep_angle_ht,A_h,S_e):
    W_dg = W_dg * kg_to_lbs
    L_t = L_t * m_to_ft
    F_w = F_w * m_to_ft
    B_h = B_h * m_to_ft
    S_e = S_e * m2_to_ft2
    S_ht = S_ht * m2_to_ft2
    K_y = 0.3 * L_t
    W_ht = 0.0379 * K_uht * (1+F_w/B_h)**(-0.25) * W_dg**(0.639) * N_z**(0.1) * S_ht**(0.75) * L_t**(-1) * K_y**(0.704) * (cos(Sweep_angle_ht))**(-1) * A_h**(0.166) * (1+S_e/S_ht)**(0.1)
    return W_ht/kg_to_lbs

def get_weight_Engine_Controls(N_en,L_ec):
    L_ec = L_ec * m_to_ft
    W_ec = 5*N_en + 0.8 * L_ec
    return W_ec/kg_to_lbs

def get_weight_furnishing(N_c,W_c,S_f):
    S_f = S_f * m2_to_ft2
    W_c = W_c * kg_to_lbs
    W_furn = 0.0577 * N_c**(0.1) * W_c**(0.393) * S_f**(0.75)
    return W_furn/kg_to_lbs

def get_weight_handling_gear(W_dg):
    W_dg = W_dg * kg_to_lbs
    return 3*10**(-4)* W_dg / kg_to_lbs

def get_weight_instruments(K_r,K_tp,N_c,N_en,L_fuselage_whole,B_w):
    L_fuselage_whole = L_fuselage_whole * m_to_ft
    B_w = B_w * m_to_ft
    W_instru = 4.509 * K_r * K_tp * N_c**(0.541) * N_en*(L_fuselage_whole+B_w)**(0.5)
    return W_instru/kg_to_lbs

def get_weight_airconditioning(N_p,V_pr,W_uav):
    W_uav = W_uav*kg_to_lbs
    V_pr = V_pr * m3_to_ft3
    W_airco = (62.36 * N_p ** 0.25 * (V_pr / 1000) ** 0.604 * W_uav ** 0.10)
    return W_airco/kg_to_lbs

def get_weight_electrical(R_kva,L_a,N_gen):
    L_a = L_a * m_to_ft
    W_electrical = (7.291 * R_kva **(0.782) * L_a ** 0.346 * N_gen ** 0.1)
    return W_electrical/kg_to_lbs
#def get_weight_nacelle_group(K_ng,N_Lt,N_w,W_ec,N_en,S_n):

def get_weight_apuinstalled(W_APU_uninstalled):
    W_APU_uninstalled = W_APU_uninstalled * kg_to_lbs
    W_apuins = 2.2 * W_APU_uninstalled
    return W_apuins/kg_to_lbs

def get_weight_starter(N_en,W_en):
    W_en = W_en*kg_to_lbs
    W_starter = 49.19 * (N_en*W_en/1000)**(0.541)
    return W_starter/kg_to_lbs

def get_weight_hydraulics(N_f,L_fuselage_whole,B_w):
    L_fuselage_whole = L_fuselage_whole * m_to_ft
    B_w = B_w * m_to_ft
    W_hydraulics = 0.2673 * N_f *(L_fuselage_whole + B_w)**(0.937)
    return W_hydraulics/kg_to_lbs

# -----------------Inputs-----------------------
L_t =  3     #Tail length [m]
W_dg = 8618             #Design gross weight [kg]
N_z = 4.3913             #Ultimate load factor (1.5* limit load factor)
B_w = 20.12                #Wing span [m]

#Wing inputs:
S_w = 45             #Wing surface [m^2]
A = 9                #Aspect ratio
t_c_root =   0.16       #Thickness to chord ratio at root [m]
Lambda =  0.5          #Wing taper ratio
Sweep_angle = 0.0      #Sweep angle at 25% MAC [rad]
S_csw = 2.54           #Control surface area [m^2]

#Avionics inputs
W_uav = 600/kg_to_lbs   #Uninstalled avionics weight

#Landing gear inputs
K_mp = 1         #1.126 for kneeling gear;= 1.0 otherwise
K_np = 1         #1.15 for kneeling gear;= 1.0 otherwise
N_mw =  2        #Number of main wheels
N_nw =  1        #number of nose wheels
N_mss = 2        #Number of main gear shock struts   (No clue how much)
V_stall = 35     #Stall velocity
L_n =   2       #nose gear length [m]
L_m =    2       #Length of main landing gear [m]
N_l = 4.3913         #Ultimate landing load factor       (Used N_Z for now)
W_l = 8618          #Landing design gross weight [kg]    (took 5% of MTOW for now)

#Fuselage inputs
K_door = 1     #1.0 if no cargo door; = 1.06 if one side cargo door; = 1.12 if two side cargo doors; = 1.12 if aft clamshell door; = 1.25 if two side cargo doors and aft clamshell door
K_lg = 1.12    # 1.12 if fuselage-mounted main landing gear;= 1.0 otherwise
L_fus =   8      #Fuselage length excluding tail cap en radome
S_f =  780 / m2_to_ft2
L_over_D = 9.2   #Lift over drag

#Vertical tail inputs
H_t = 1.5/4     #Horizontal tail height above fuselage [m]
H_v = 1     # Vertical tail height above fuselage [m]
S_vt = 6.7 * (45/35.18)    #Surface area vertical tail [m^2]
Sweep_angle_vt =   (40/180)*pi     #Sweep angle vertical tail [rad]
A_v =  1.58    #Aspect ratio vertical tail

#Inputs Horizontal tail
K_uht = 1           #1.143 for unit (all-moving) horizontal tail; = 1.0 otherwise
F_w =  1.1             #fuselage width at horizontal tail intersection, [m]   (Used the finger method of Timo)
B_h = 6.75 * (20.12/19.48)         # Horizontal tail span [m]
S_ht = 9.88 * (45 / 35.18)          #horizontal tail surface [m^2]
A_h = 6.73 * (9 / 11.45)        #Aspect ratio of horizontal tail
Sweep_angle_ht = (8/180)*pi     #Horizontal tail sweep [rad]
S_e = 1            #Elevator area [m^2]

#Inputs engine control
N_en = 2        #Is number of engines
L_ec = 4.7         #length from engine front to cockpit-total if multiengine, [m]

#Inputs furnishing
N_c = 2         #Number of crew members
W_c = 1814      #Max cargo (payload?)

#Instrument inputs
K_r = 1         #1.133 if reciprocating engine; = 1.0 otherwise
K_tp = 1        # = 1
L_fuselage_whole = 11.028

#Airconditioning inputs
N_p = 21        #Nr of crew + number of passengers
V_pr = 632/m3_to_ft3

#Electrical inputs
R_kva=	60     # system electrical rating, [kv · A]
L_a= 6.65       #electrical routing distance, generators to avionics to cockpit [m]
N_gen = 2       #number of generators (typically = N_en)

#APU inputs
W_APU_uninstalled = 36      #Weight uninstalled auxiliary power unit

#Starter inputs
W_en =  250        #Engine weight [kg]

#Inputs hydraulics
N_f = 5     #Number of control surfaces (2x aileron, 2x elevator, 1x rudder)




weight_wing = get_weight_wing(W_dg,N_z,S_w,A,t_c_root,Lambda,Sweep_angle,S_csw)
weight_avionics = get_weight_avionics(W_uav)
weight_landing_gear = get_weight_landing_gear(K_mp,K_np,N_l,W_l,L_m,L_n,N_mw,N_nw,N_mss,V_stall)
weight_fuselage = get_weight_fuselage(K_door,K_lg,W_dg,N_z,L_fus,S_f,L_over_D,B_w,Lambda)
weight_vertical_tail = get_weight_vert_tail(H_t,H_v,W_dg,N_z,L_t,S_vt,Sweep_angle_vt,A_v,t_c_root)
weight_horizontal_tail = get_weight_hor_tail(K_uht,F_w,B_h,W_dg,N_z,S_ht,L_t,Sweep_angle_ht,A_h,S_e)
weight_engine_control = get_weight_Engine_Controls(N_en,L_ec)
weight_furnishing = get_weight_furnishing(N_c,W_c,S_f)
weight_handling_gear = get_weight_handling_gear(W_dg)
weight_instruments = get_weight_instruments(K_r,K_tp,N_c,N_en,L_fuselage_whole,B_w)
weight_airconditioning = get_weight_airconditioning(N_p,V_pr,W_uav)
weight_electrical = get_weight_electrical(R_kva,L_a,N_gen)
weight_APUins = get_weight_apuinstalled(W_APU_uninstalled)
weight_starter = get_weight_starter(N_en,W_en)
weight_hydraulics = get_weight_hydraulics(N_f,L_fuselage_whole,B_w)

print("Weight wing =", weight_wing, "   % of MTOW:", 100*weight_wing/W_dg)
print("Weight avionics =", weight_avionics, "   % of MTOW:", 100*weight_avionics/W_dg)
print("Weight Landing gear =", weight_landing_gear, "   % of MTOW:", 100*weight_landing_gear/W_dg)
print("Weight fuselage = ", weight_fuselage, "   % of MTOW:", 100*weight_fuselage/W_dg)
print("Weight vertical tail = ", weight_vertical_tail, ",   % of MTOW:", 100*weight_vertical_tail/W_dg)
print("Weight horizontal tail = ", weight_horizontal_tail, ",   % of MTOW:", 100*weight_horizontal_tail/W_dg)
print("Weight engine control =", weight_engine_control, ",   % of MTOW:", 100*weight_engine_control/W_dg)
print("Weight furnishing =", weight_furnishing, ",   % of MTOW:", 100*weight_furnishing/W_dg)
print("Weight handling gear =", weight_handling_gear, ",   % of MTOW:", 100*weight_handling_gear/W_dg)
print("Weight instruments =", weight_instruments, ",   % of MTOW:", 100*weight_instruments/W_dg)
print("Weight airconditioning =", weight_airconditioning, ",   % of MTOW:", 100*weight_airconditioning/W_dg)
print("Weight electrical =", weight_electrical, ",   % of MTOW:", 100*weight_electrical/W_dg)
print("Weight installed Auxiliary power unit=", weight_APUins, ",   % of MTOW:", 100*weight_APUins/W_dg)
print("Weight starter (pneumatic) =", weight_starter, ",   % of MTOW:", 100*weight_starter/W_dg)
print("Weight hydraulics =", weight_hydraulics, ",   % of MTOW:", 100*weight_hydraulics/W_dg)
print("\n Total weight of subsystems:", weight_wing + weight_avionics + weight_landing_gear + weight_fuselage + weight_vertical_tail + weight_horizontal_tail + weight_engine_control + weight_furnishing + weight_handling_gear + weight_instruments + weight_airconditioning + weight_electrical + weight_APUins + weight_starter + weight_hydraulics)