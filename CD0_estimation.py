from math import *

ft = 3.2808399
ft2 = 10.7639104
ft3 = 35.3146667
lbs = 2.20462262
inch = 39.3700787
lbft2 = 23.73025

def get_S_w(S_w_w,S_w_f,S_w_t,S_w_n):
    S_w = S_w_w + S_w_f + S_w_t + S_w_n
    return S_w

def get_S_w_w(S, b, f_d, lambda_t, tc_max):
    lambda_b = 1 + f_d/b(lamda_t - 1)
    S_e_w = 1 - (f_d/b*(1+lambda_b))/(lambda_t+1)
    pc = 2*(1+tc_max**2)
    S_wing = S*pc*S_e_w
    return S_wing

def get_S_w_f(l_f, l_tc, l_nc, d)
    S_fuselage = pi*d*l_f(1-1/3*(l_nc/l_f)-1/2*(l_tc/l_f))
    return S_fuselage

def get_S_w_t(S)
    if S*ft2 < 3000:
        S_tail = 0.88*S
    else:
        S_tail = 2.5*S^0.85
    return S_tail

def get_S_w_n(d_n,l_n)
    S_nacelle = 2*pi*d_n*l_n
    return S_nacelle



####### inputs ######

# Wing/tail inputs #
b = 20*ft                                   #Wing span [m]
yb = 1*ft                                     #
S = 40*ft2                                      #Wing area [m]
tc_max = 0.17
lambda_t = 0.5

V = 90                                      #cruise speed [m/s]
rho = 0.9                                   #density
mu = 14.16*10**-6                           #kinematic viscosity [m^2/s] @10 Degrees Celsius

# fuselage inputs #
l_f = 12.24
l_tc = 4.78
l_nc = 3
f_d = 2.2

# nacelle inputs #
d_n = 0.3
l_n = 1.8

#--------------------------------------------------#
S_w = 1/ft2*get_S_w(S_w_w,S_w_f,S_w_t,S_w_n)      #Total wetted area in clean configuration
S_w_w = 1/ft2*get_S_w_w(S, b, f_d, lambda_t, tc_max)                         #Wing wetted area
S_w_f = 1/ft2*get_S_w_f(l_f, l_tc, l_nc, d)                                   #Fuselage wetted area
S_w_n = 1/ft2*get_S_w_n(d_n,l_n)                                   #Nacelle wetted area
S_w_t = 1/ft2*get_S_w_t(S)                        #Tail wetted area

R_e = (rho*V/mu)*(S_w/b)
Cd0 = 0.00258+0.00102*e**(-6.28*10**-9*R_e)+0.00295*e**(-2.01*10**-8*R_e)*(S_w/b)/(S/b)

print("CD_0 ="Cd0)