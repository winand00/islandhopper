from math import *

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
    S_e_w_2

def get_S_w_w()


S_wing =

S_w = 1     #Total wetted area in clean configuration
S_w_w = 1   #Wing wetted area
S_w_f = 1   #Fuselage wetted area
S_w_n = 1   #Nacelle wetted area

#Wing wetted area
