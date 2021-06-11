from math import *

m_to_ft = 3.2808399
m2_to_ft2 = 10.7639104
m3_to_ft3 = 35.3146667
kg_to_lbs = 2.20462262
m_to_in = 39.3700787
kgm2_to_lbft2 = 23.73025

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



S_w = 1     #Total wetted area in clean configuration
S_w_w = 1   #Wing wetted area
S_w_f = 1   #Fuselage wetted area
S_w_n = 1   #Nacelle wetted area

#Wing wetted area
b = 20.12       #Wing span [m]
c_r = 1.5         #Chord at root [m]
c_t = 2.25         #Chord at tip [m]
y_c = 1.92/2         #Distance from mid of wing to where wing meets fuselage [m] (0.5 times diameter)
