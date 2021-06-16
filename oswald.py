# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 15:04:42 2021

@author: Anthony Cummins
"""

#cd0 = 0.041
# M = 0.274
import numpy as np



flambda = lambda x: 0.0524*x**4 -0.15*x**3 + 0.1659*x**2 - 0.0706*x + 0.0119

dlambda = -0.357 + 0.45

def e_theo(lmbda=0.5, dlambda=dlambda, A = 8.889):
    return 1/(1+flambda(lmbda - dlambda)*A)


def e1(e_theo, k_ef, k_edo, k_em):
    return e_theo*k_ef*k_edo*k_em

def k_ef(df=2.255, b=20):
    return 1- 2*(df/b)**2


def e2(C_d0=0.0252, A=8.889):
    K = 0.38
    P = K*C_d0
    Q = 1/e_theo()*k_ef()
    k_em =1
    
    return k_em/ (Q+P*np.pi*A)

print(e2())
