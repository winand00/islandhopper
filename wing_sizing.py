# -*- coding: utf-8 -*-
"""
Created on Fri May 21 10:41:19 2021

@author: Anthony Cummins

AILERON EFFECTIVENESS TAU MUST BE MANUALLY SET FOR ITERATION
"""


import numpy as np
import scipy as sp
import scipy.integrate
import math


#design Cl
m = 8618
g = 9.81
S = 45
W = m*g
V_cruise = 90
V_ldg = 55
rho_cruise = 0.9046
rho_sl = 1.225
q_cruise = 0.5*rho_cruise*V_cruise*V_cruise
A = 9
tc = 0.16 #thickness-to-chord ratio, function of Mach no.
b = np.sqrt(S*A)
c = S/b
mu= 1.84*10**-5

Re_cruise = rho_cruise*V_cruise*c/mu
Re_ldg = rho_sl*V_ldg*c/mu

#size airfoil based on design lift 
#Raymer p47 for t/c ratio based on Mach

Cl_des = 1.1*W/q_cruise/S #design lift coefficient

Cl_max_airfoil = 1.8 #maximum lift coeff of airfoil, from xfoil

Cl_max_wing = 0.9*Cl_max_airfoil #maximum lift coeff of wing

Cl_max_ldg = 2.4 #landing cl

dCL_max =Cl_max_ldg - Cl_max_wing #increment in lift coeff

#Aileron
tau = 0.3 #aileron effectiveness, function of aileron-chord-to-total-chord ratio
lamb = 0.5 #taper ratio
clalpha = 0.127  # airfoil lift curve slope, from xfoil
#b1 and b2 middle of plane is 0
cd0 = 0.006 #drag coefficient of AIRFOIL
b1 = 0.7*b/2 #start of aileron
b2 = 0.9*b/2 #end of aileron

y = np.arange(0,b/2+0.01,0.01)





def flap_calculator(dCl_max_TE, dCL_max= dCL_max, LE = False): 
    """
    

    Parameters
    ----------
    dCl_max_TE : This is the increment in Cl_max of the TE FLAP. Select yourself.
        
    dCL_max : This is the required increase in CL of the ENTIRE WING
    The default is dCL_max, which has been calculated above
    
    LE : States whether LE devices are present or not. The default is False. 
    If LE devices are present, LE is True. The function will prompt for the 
    increase in Cl of the Leading Edge AND the Trailing Edge devices.

    Returns
    -------
    SwfS_LE : Flapped surface fraction of the leading edge
       
    SwfS_TE : Flapped surface fraction of the trailing edge

    """
    if LE == False:
        dCL_max_TE = dCL_max
        SwfS_TE = dCL_max_TE/0.9/dCl_max_TE
        SwfS_LE = 0
        
    elif LE == True:
        dCL_max_LE = float(input("Give the dCL_max of leading edge device: "))
        dCl_max_LE = float(input("Give the dCl_max of leading edge device: "))
        SwfS_LE = dCL_max_LE/0.9/dCl_max_LE
        dCL_max_TE = dCL_max - dCL_max_LE
        SwfS_TE = dCL_max_TE/0.9/dCl_max_TE
        
    return SwfS_LE, SwfS_TE


def chord(y):
    return (2*S/((1+lamb)*b))*(1-(1-lamb)*abs(2*y)/b)

def chordy(y):
    return (2*S/((1+lamb)*b))*(1-(1-lamb)*abs(2*y)/b)*y

def chordy2(y):
    return (2*S/((1+lamb)*b))*(1-(1-lamb)*abs(2*y)/b)*y*y
    


def roll_rate(clalpha=clalpha, tau=tau, S=S, b=b, b1=b1, b2=b2):

    cl_da = 2*clalpha*tau/(S*b)*sp.integrate.quad(chordy,b1,b2)[0] #due to aileron def.
    cl_p = -4*(clalpha+cd0)/(S*b*b)*sp.integrate.quad(chordy2,b1,b2)[0] #roll damping
    
    da = math.radians(25) #aileron deflection, radians
    P = math.degrees(-cl_da/cl_p*da*(2*V_ldg/b))
    
    return P


flapped_area = flap_calculator(1.3)
print("flapped area fraction = ",flapped_area)
rollrate = roll_rate()
print("roll rate = ", rollrate)
 




