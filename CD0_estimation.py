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

def get_S(b,c_r,c_t,y_c,c_c):
    Lambda = c_t/c_r
    Lambda_c = c_c/c_r
    eta_c = y_c * (b/2)
    S = 0.5 * b * c_r * (Lambda*(1-eta_c) + eta_c + Lambda_c)
    return S

def get_S_e_w(S,c_r,d,):
    S_e_w = S - c_r * d
    S_e_w_2 =

def get_S_w_w()

def get_S_w_f(l_f, l_tc, l_nc, d)
    S_fuselage = pi*d*l_f(1-1/3*(l_nc/l_f)-1/2*(l_tc/l_f))
    return S_fuselage

def get_S_w_t(S)
    if S*10.7639 < 3000:



S_w = 1     #Total wetted area in clean configuration
S_w_w = 1   #Wing wetted area
S_w_f = 1   #Fuselage wetted area
S_w_n = 1   #Nacelle wetted area

#Wing wetted area
b = 20.12       #Wing span [m]
c_r = 1.5         #Chord at root [m]
c_t = 2.25         #Chord at tip [m]
y_c = 1.92/2         #Distance from mid of wing to where wing meets fuselage [m] (0.5 times diameter)
