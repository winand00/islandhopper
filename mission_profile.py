import numpy as np
from ISA import script
import matplotlib.pyplot as plt

P_low = 25000 #watt
Pmax = 1150000
V_climb = 60 #TBD
tclimb = 10*60

#constants
E_hydrogen = 120000000 #J/kg
efficiency_constant = 1.23
eff_engine = 0.93
eff_propeller = 0.85
eff_pmad = 0.989
g = 9.81
e = 0.7
A = 10
m = 8618
W = m*g
rho_climb = 1.05
S = 40
Cd0 = 0.03
Hcruise = 3048

CL_takeoff = 2.4
rho_takeoff = 1.225
s_takeoff = 750

rho_cruise = script(Hcruise)
V_cruise = 90
distance = 555000

rho_loiter = 1
t_loiter = 45*60 #min

V_descent = 90 #TBD
descent_angle = -3 #deg
rho_descent = 1


def CD(Cl):
    return Cd0 + Cl**2/(np.pi*A*e)

def CL(rho,V):
    return W/(0.5*rho*V**2*S)

def Velocity(rho,Cl):
    return np.sqrt(2*W/(rho*S*Cl))

def Energy(t,P, eff_fuelcell):
    return t*P/(eff_propeller*eff_engine*eff_fuelcell*eff_pmad)

def efficiency(A):
    V = -0.047 * np.log(A) + 0.9782
    return V/efficiency_constant

#Take-off
#V_takeoff = Velocity(rho_takeoff,CL_takeoff)
#CD_takeoff = CD(CL_takeoff)
t_takeoff = 42
t_climbout = 132

#v = 3
#t = 0
#dt = 0.01
#fric = 0.03
#s = 0

#while v < V_takeoff:
 #   T = Pmax/v
  #  D = CD(CL_takeoff)*0.5*rho_takeoff*v**2*S
   # Ffric = fric*W
    #a = (T-D-Ffric)/(m)
  #  v = v + a*dt
  #  ds = v*dt
  #  s = s +ds
  #  t = t+dt
#print("runway =",s)


P_takeoff = Pmax+P_low
P_climbout = 0.85*P_takeoff
RC_climbout = 4
Hclimbout = t_climbout*RC_climbout

#climb
Cl_climb = CL(rho_climb,V_climb)
Cd_climb = CD(Cl_climb)
Pr = Cd_climb/Cl_climb *W*V_climb

#RC = (Hcruise-Hclimbout)/tclimb
RC = V_climb*(np.sin(np.radians(1.37)))

P_climb = RC*W + Pr

#climb one engine inoperative case 1
rho_climbOEI = 1.225
V_climbOEI = 1.2 * Velocity(rho_climbOEI, 1.8)
print("Speed", V_climbOEI)
CL_climbOEI = CL(rho_climbOEI, V_climbOEI)
print("CL", CL_climbOEI)
Cd_climbOEI = CD(CL_climbOEI)
Pr = Cd_climbOEI/Cl_climb *W*V_climbOEI
print("Power required", Pr)
RC_OEI = V_climbOEI*(np.sin(np.radians(1.14576)))
P_climbOEI = RC_OEI*W + Pr
P_engineOEI = P_climbOEI / eff_propeller
print("One engine inoperative", P_climbOEI, P_engineOEI)

#climb one engine inoperative case 1
rho_climbOEI2 = 1.17215
V_climbOEI2 = 1.2 * Velocity(rho_climbOEI2, 1.8)
CL_climbOEI2 = CL(rho_climbOEI2, V_climbOEI2)
Cd_climbOEI2 = CD(CL_climbOEI2)
PrOEI2 = Cd_climbOEI2/CL_climbOEI2 *W*V_climbOEI2
RC_OEI2 = V_climbOEI*(np.sin(np.radians(0.6875)))
P_climbOEI2 = RC_OEI2*W + PrOEI2
P_engineOEI2 = P_climbOEI2 / eff_propeller
print("One engine inoperative case 2", P_climbOEI2, P_engineOEI2)


#cruise
Cl_cruise = CL(rho_cruise, V_cruise)
Cd_cruise = CD(Cl_cruise)
tcruise = distance/V_cruise

Pr_cruise = Cd_cruise/Cl_cruise *W *V_cruise

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
t_descent = -Hcruise/(V_descent*np.sin(np.radians(descent_angle)))

P_descent = T*V_descent

#sizing

#battery
power_density = 2000 #W/kg
energy_density = 360 *3600#Wh/kg


inputfactor = 1/(eff_propeller*eff_engine*eff_pmad)

P_fuelcell = Pr_cruise*inputfactor +P_low
m_fuelcell = 111.86*np.log(P_fuelcell/1000) - 261.95 + P_fuelcell/3000
Pbat = P_takeoff*inputfactor-P_fuelcell
mbat = Pbat/power_density
Ebat = mbat * energy_density

print("Battery until what phase?: climbout/climb")
batphase = input()

Ebat_needed = (P_takeoff*inputfactor-P_fuelcell)*t_takeoff + (P_climbout*inputfactor-P_fuelcell)*t_climbout + (Pmax*inputfactor*0.1*10*60+P_low)
if batphase == "climb":
    Ebat_needed= Ebat_needed + (P_climb*inputfactor-P_fuelcell)*tclimb

if Ebat_needed>Ebat:
    mbat = Ebat_needed/(energy_density)
    Pbat = mbat * power_density
    print("Energy limiting")
else:
    print("power limiting")

print("Fuel cell power = ",P_fuelcell/1000, "kW")
print("Fuel cell mass",m_fuelcell )
print("Battery power = ",Pbat/1000, "kW")
print("Battery energy = ", mbat*energy_density)
print("Battery mass = ", mbat)

#Hydrogen calculation
flight_time = t_takeoff+tclimb+tcruise+t_descent+t_loiter
E_takeoff = Energy(t_takeoff,P_takeoff,0.53) + Energy(t_climbout,P_climbout,0.55)
E_climb = Energy(tclimb,P_climb,0.57)
E_cruise = Energy(tcruise,Pr_cruise,0.585)
E_descent = Energy(t_descent,P_descent,0.63)
E_loiter = Energy(t_loiter,P_loiter,0.61)
E_taxi = Energy(15*60,Pmax*0.1,0.67)

Etotal = E_takeoff +E_climb+E_cruise+E_descent+E_taxi+E_loiter
Hydrogen = Etotal/E_hydrogen
print(Hydrogen)
print('fight time = ', flight_time/3600)

#radiator
Reynolds = 15000000
T_diff = 50
Pratzl = 0.71
#Nusselt =
print(E_taxi/3600000,E_takeoff/3600000,E_climb/3600000,E_cruise/3600000,E_descent/3600000, Etotal/3600000)

