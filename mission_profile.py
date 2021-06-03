import numpy as np
from ISA import script
import matplotlib.pyplot as plt

g = 9.81
V_climb = 60 #TBD
e = 0.8
A = 9
m = 8618
W = m*g
rho_climb = 1.05
S = 45
Cd0 = 0.03
tclimb = 10*60
Hcruise = 3048

CL_takeoff = 2
rho_takeoff = 1.225
s_takeoff = 750

rho_cruise = script(Hcruise)
V_cruise = 90
distance = 370400

rho_loiter = 1
t_loiter = 45 #min

V_descent = 90 #TBD
descent_angle = -3 #deg
rho_descent = 1


#battery
power_density = 5000 #W/kg
energy_density = 110 #Wh/kg




def CD(Cl):
    return Cd0 + Cl**2/(np.pi*A*e)

def CL(rho,V):
    return W/(0.5*rho*V**2*S)

def Velocity(rho,Cl):
    return np.sqrt(2*W/(rho*S*Cl))

#Take-off
V_takeoff = Velocity(rho_takeoff,CL_takeoff)
CD_takeoff = CD(CL_takeoff)

v = 3
t = 0
dt = 0.01
Pmax = 1300000
fric = 0.03
s = 0

while v < V_takeoff:
    T = Pmax/v
    D = CD(CL_takeoff)*0.5*rho_takeoff*v**2*S
    Ffric = fric*W
    a = (T-D-Ffric)/(m)
    v = v + a*dt
    ds = v*dt
    s = s +ds
    t = t+dt
#print("runway =",s)

P_takeoff = 1300000
P_climbout = 0.85*P_takeoff


#climb
Cl_climb = CL(rho_climb,V_climb)
Cd_climb = CD(Cl_climb)
Pr = Cd_climb/Cl_climb *W*V_climb

RC = 3000/tclimb

P_climb = RC*W + Pr
#print("climb power =", P_climb)
#print("RC =", RC)
#cruise
Cl_cruise = CL(rho_cruise, V_cruise)
Cd_cruise = CD(Cl_cruise)
Pr_cruise = Cd_cruise/Cl_cruise *W *V_cruise

tcruise = distance/V_cruise
energyused = tcruise*Pr_cruise/3600000

#loiter
Cl_loiter = np.sqrt(3*Cd0*np.pi*A*e)
V_loiter = Velocity(rho_loiter,Cl_loiter)
Cd_loiter = CD(Cl_loiter)

D = Cd_loiter/Cl_loiter * W
P_loiter = D*V_loiter

#decent
Cl_descent = CL(rho_descent,V_descent)
Cd_descent = CD(Cl_descent)

D = Cd_descent*0.5*rho_descent*V_descent**2*S
T = D + W*np.sin(np.radians(descent_angle))

#P_descent =

#sizing
t_takeoff = 42
t_climbout = 132
P_fuelcell = 500000
Pbat = P_takeoff-P_fuelcell
mbat = Pbat/power_density
Ebat = mbat *energy_density

Ebat_needed = ((P_takeoff-P_fuelcell)*t_takeoff + t_climbout*(P_climbout-P_fuelcell)+(P_climb-P_fuelcell)*tclimb)/3600

if Ebat_needed>Ebat:
    mbat = Ebat_needed/energy_density
    Pbat = mbat *power_density
    print("Energy limiting")
else:
    print("power limiting")

print("Fuel cell power = ",P_fuelcell/1000, "kW")
print("Battery power = ",Pbat/1000, "kW")
print("Battery mass = ", mbat)



