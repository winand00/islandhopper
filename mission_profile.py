import numpy as np

g = 9.81
V_climb = 100
e = 0.8
A = 9
W = 8600*g
rho_climb = 1.05
S = 45
Cd0 = 0.03
tclimb = 10*60
Hcruise = 3000

def CD(Cl):
    return Cd0 + Cl**2/(np.pi*A*e)

def CL(rho,V):
    return W/(0.5*rho*V**2*S)

#Take-off



#climb
Cl_climb = CL(rho_climb,V_climb)
Cd_climb = CD(Cl_climb)
Pr = Cd_climb/Cl_climb *W*V_climb

RC = 3000/tclimb

Pa = RC*W +Pr

print(Pa)

#cruise



#loiter


#decent


