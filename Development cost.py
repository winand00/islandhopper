from math import *


W_ampr = 4000 * 2.20462262 #[lbs]
V_max = 175 #[kts]
N_rdte = 2 #[-]
N_st = 1
CEF = 4 #[-]
F_diff = 2 #[-]
F_cad = 0.8 #[-]


C_dstr = 0.008325 * (W_ampr ** 0.873) * (V_max**1.890) * (N_rdte**0.346) * CEF * F_diff
print(C_dstr)

MHR_aedr = 0.0396 * (W_ampr ** 0.791) * (V_max**1.526) * (N_rdte**0.183) * F_diff * F_cad
print(MHR_aedr)

C_ftar = 7.8 * 10 ** 6 * 0.75 #nog niet zeker
print(C_ftar)

C_ftor = 0.001244 * (W_ampr ** 1.160) * (V_max ** 1.371) * (N_rdte - N_st)**1.281 * CEF * F_diff
print(C_ftor)

C_tsfr = 0 #add 20 \% for extra testing facilities to final cost

C_pror = 0 #add 10 \% profit to final cost

C_finr = 0 #add 10 \% interest to final cost

print('cost: ', C_dstr + MHR_aedr*150 + C_ftar + C_ftor)


