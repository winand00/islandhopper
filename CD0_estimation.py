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

def get_S_w_w(S, b, y_b, lambda_t, tc_max):
    lambda_b = 1. + y_b/b*(lambda_t - 1.)
    S_e_w = 1. - (y_b/b*(1.+lambda_b))/(lambda_t+1.)
    pc = 2*(1+tc_max**2)
    S_wing = S*pc*S_e_w
    return S_wing

def get_S_w_f(l_f, l_tc, l_nc, f_d):
    S_fuselage = pi*f_d*l_f*(1-1/3*(l_nc/l_f)-1/2*(l_tc/l_f))
    return S_fuselage

def get_S_w_t(S):
    if S < 3000:
        S_tail = 0.88*S
    else:
        S_tail = 2.5*S**0.85
    return S_tail

def get_S_w_n(d_n,l_n):
    S_nacelle = 2*pi*d_n*l_n
    return S_nacelle



####### inputs ######

# Wing/tail inputs #
b = 20.*ft                                   #Wing span [m]
y_b = 0.5*ft                                     #
S = 40.*ft2                                      #Wing area [m]
tc_max = 0.17
lambda_t = 0.5

V = 90.                                     #cruise speed [m/s]
rho = 0.9                                 #density
mu = 14.16*(10E-6)                           #kinematic viscosity [m^2/s] @10 Degrees Celsius

# fuselage inputs #
l_f = 12.24*ft
l_tc = 4.78*ft
l_nc = 3*ft
f_d = 2.2*ft

# nacelle inputs #
d_n = 0.3*ft
l_n = 1.8*ft

#--------------------------------------------------#

S_w_w = get_S_w_w(S, b, y_b, lambda_t, tc_max)      #Wing wetted area
S_w_f = get_S_w_f(l_f, l_tc, l_nc, f_d)             #Fuselage wetted area
S_w_n = get_S_w_n(d_n,l_n)                          #Nacelle wetted area
S_w_t = get_S_w_t(S)                                #Tail wetted area
S_w = S_w_w + S_w_f + S_w_t + S_w_n                 #Total wetted area in clean configuration
Swb = 10.7*(S/b)**0.75

R_e = (rho*V/mu)*(Swb/ft)
CFe = 0.00258+0.00102*e**(-6.28E-9*R_e)+0.00295*e**(-2.01E-8*R_e)
Cd0 = CFe*10.7*(S/b)**(0.75-1)

print("S_w_w =", S_w_w/ft2)
print("S_w_f =", S_w_f/ft2)
print("S_w_n =", S_w_n/ft2)
print("S_w_t =", S_w_t/ft2)
print("CD_0 =", Cd0)
print(Swb)
