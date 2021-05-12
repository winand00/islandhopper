import numpy as np
from math import *

def prop_diameter(N,dY, dy, b):
    return dY*b/(N*(1+dy))

def normalized_D_p(N,dY, dy,A,W_S):
    return dY**2/(N**2*(1+dy)**2) * A/(W_S)

def thrust_coef(T_W,N,rho, V, D_normalized, chi):
    return chi*T_W/(rho*V**2*D_normalized*N)

def a_p(T_c):
    return 0.5*(sqrt(1+8*T_c/pi)-1)

def contraction_ratio(a_p, xp_rp):
    return sqrt((1+a_p)/(1+a_p*(1+xp_rp)/sqrt(xp_rp**2+1)))

def a_w(a_p, contraction):
    return (a_p+1)/contraction**2 - 1

def alpha_w(M,A,sweep, CL):
    return CL/(2*pi*A)*(2+sqrt(A**2*(1-M**2)*(1+tan(sweep)**2/(1-M**2))+4))

def dcl(a_w,alpha_w, alpha_p, beta):
    return 2*pi*((sin(alpha_w)-a_w*beta*sin(alpha_p-alpha_w))*sqrt((a_w*beta)**2 + 2*a_w*beta*cos(alpha_p) + 1) - sin(alpha_w))





def DCL(dY,M,A,sweep,CL,alpha_p,beta,xp_rp, T_W,rho,V,N,dy,W_S,chi):
    return dcl(a_w(a_p(thrust_coef(T_W,N,rho, V,normalized_D_p(N,dY, dy,A,W_S),chi)),contraction_ratio(a_p(thrust_coef(T_W,N,rho, V,normalized_D_p(N,dY, dy,A,W_S), chi)),xp_rp)), alpha_w(M,A,sweep,CL), alpha_p, beta)*dY

def DCD0(dY, xp_rp,cf, T_W,rho,V,N,dy,W_S,chi,A):
    return dY*a_w(a_p(thrust_coef(T_W,N,rho, V,normalized_D_p(N,dY, dy,A,W_S),chi)),contraction_ratio(a_p(thrust_coef(T_W,N,rho, V, normalized_D_p(N,dY, dy,A,W_S), chi)),xp_rp))**2*cf

dY = 0.7
V = 90 #m/s
a = 328.4
M = V/a
A = 15
sweep = 0
CL = 0.5
alpha_p = 0.
beta = 0.8
xp_rp = 0.7
T_W = 0.645
rho = 0.904637
cf = 0.009
e = 0.8
cd0 = 0.03

N = 12
dy = 0.05
W_S = 2842
chi = 1
b = 20

dcl = DCL(dY,M,A,sweep,CL,alpha_p,beta,xp_rp, T_W,rho,V,N,dy,W_S,chi)
dcd0 = DCD0(dY, xp_rp,cf, T_W,rho,V,N,dy,W_S,chi,A)
dcdi = (dcl**2 + 2*CL*dcl)/(pi*A*e)


oldcl = CL
oldcd = cd0 + CL**2/(pi*A*e)

clopt = sqrt(1/3 *cd0*A*e)

newcl = CL+dcl
newcd = oldcd+dcdi+dcd0
newcd0 = cd0+dcd0


print(oldcl,newcl)
print(oldcd,newcd)
print(oldcl/oldcd,newcl/newcd)

print(prop_diameter(N,dY, dy, b))
print(newcd0)