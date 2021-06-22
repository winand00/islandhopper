# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 12:32:39 2021

@author: Anthony Cummins
"""

import numpy as np
from matplotlib import pyplot as plt

# =============================================================================
# 
# def r(x, x_0, x_1, r_0, r_1):
#     b = (r_1-r_0)/(x_1-x_0)
#     r = b*(x-x_0)+ r_0
#     return r
# =============================================================================
    
    

def A(r):   
    A = np.pi*r**2
    return A
    

def V(A_0, A_1, V_0):
    
    V_1 = (A_0/A_1)*V_0
    return V_1

def T(T_0, V_0, A_0, A_1, C_p, beta):
    
    gamma = beta**2
    
    zeta = (1-gamma*(A_0/A_1)**2)/2/C_p
    
    T_1 = T_0 + zeta*V_0**2
    
    return T_1

def Mach(gamma, R, T, V):
    a = np.sqrt(gamma*R*T)
    M = V/a
    return M

def p (beta, rho, R, T):
    
    p = rho/beta*R*T
    
    return p

def beta(M, gamma):
    beta = (1+ (gamma-1)/2*M**2)**(-1/(gamma-1)) 
    return beta


def calcQdot(h,A_rad,T_inf, T_rad):
    T_s = (T_rad+T_inf)*0.5
    Qdot = h*A_rad*(T_inf-T_s)
    return Qdot


def calcRe(rho, V, L, mu):
    return rho*V*L/mu

def calcPr(mu, c_p, k):
    return mu*c_p/k

def calch(Re, Pr, k, L_rad):
    return 0.037*Re**0.8*Pr**0.333*k/L_rad
    

###########################################
#input parameters



#constants
R = 287
T_0 = 278 
gamma = 1.4
rho_2 = 0.9
beta_23 = 0.95
c_p = 1000
#########################################

#inlet parameters
x_2 = 0
x_3 = 1
r_2 = 0.3
r_3 = 0.6


#radiator parameters

L_rad = 1
W_rad = 1
############################################


    
#inlet calculations
A_2 = A(r_2)
A_3 = A(r_3)


V_2 = 95 
T_2 = T_0
p_2 = rho_2*R*T_0




#just before radiator section 
V_3 = V(A_2, A_3, V_2)

T_3 = T(T_2, V_2, A_2, A_3, c_p, beta_23)

#calculating static pressure 
p_3 = beta_23*rho_2*R*T_3




print("Expansion Phase Parameters: \n")
print("Flow Velocity in radiator section = ", V_3, "m/s")
print("Pressure change = ", (p_3-p_2), "Pa")
print("Temperature before radiator = ", T_3, "K")
print("")


mu = 1.48*10**-5 
k = 27*10**-3

Re = calcRe(rho_2, V_3, L_rad, mu)
Pr = calcPr(mu, c_p, k)
h = calch(Re, Pr, k, L_rad)
A_rad = L_rad*W_rad

T_inf = T_3
T_rad = 273.15+80


Qdot = calcQdot(h, A_rad, T_inf,T_rad)

n_plates = 320*10**3/Qdot


A_4 = A_3
A_5 = A_4

#temperature and velocity after radiator
T_5 = (T_rad+T_inf)/2
V_5 = np.sqrt(2*c_p*(T_5-T_3)+V_3**2) 


print("Heating phase parameters:\n")
print("Number of plates = ", int(abs(round(n_plates,0))))
print("dimensions of one plate = ",L_rad, "m x", W_rad, "m")
print("Velocity after radiator = ",V_5, "m/s \n")



r_6 = 0.5
A_6 = A(r_6)

V_6 = V(A_5, A_6, V_5)


print("Exit Velocity =", V_6, "m/s")

p_atm = 69.7*10**3

#Todo: Mach number, density relation, isentropic


    