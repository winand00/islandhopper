from math import *


K_r	= 1
K_tp=	1
N_c	=2
N_en	=2
L_f	=36.18110236
B_w	=66.01049869
R_kva=	60
L_a= 6.65 /0.0348
N_gen=	2
W_uav	=600
W_c	= 92 * 2.2 * 19 #627 * 6 # All payload
S_f	=780
N_p	=21
V_pr=	632
W_apu=	80

W_instruments = (4.509 * K_r * K_tp * N_c ** 0.541 * N_en * (L_f + B_w)**0.5) /2.2
W_electrical = (7.291 * R_kva**0.782 * L_a ** 0.346 * N_gen**0.1)/2.2
W_avionics = (1.73 * W_uav**0.983)/2.2
W_furnishings = (0.0577 * N_c**0.1 * W_c ** 0.393 * S_f ** 0.75)/2.2
W_airconditioning = (62.36 * N_p**0.25 * (V_pr/1000)**0.604 * W_uav**0.10)/2.2
W_apuins = (2.2 * W_apu)/2.2


print("W_instruments",W_instruments)
print("W_electrical",W_electrical)
print("W_avionics",W_avionics)
print("W_furnishings",W_furnishings)
print("W_airconditioning",W_airconditioning)
print("W_apuins",W_apuins)

