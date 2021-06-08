import numpy as np
from ISA import script
import matplotlib.pyplot as plt

P_low = 25000 #watt
Pmax = 1150000
V_climb = 60 #TBD
tclimb = 10*60

#constants
E_hydrogen = 120000000 #J/kg
eff_engine = 0.93
eff_propeller = 0.85
eff_pmad = 0.989
g = 9.81
e = 0.7
A = 9
m = 8618
W = m*g
rho_climb = 1.05
S = 45
Cd0 = 0.03
Hcruise = 3048

CL_takeoff = 2.4
rho_takeoff = 1.225
s_takeoff = 750

rho_cruise = script(Hcruise)
V_cruise = 90
distance = 555000

rho_loiter = 1
t_loiter = 45 #min

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
P_climbout = 0.85*P_takeoff + P_low
RC_climbout = 4
Hclimbout = t_climbout*RC_climbout

#climb
Cl_climb = CL(rho_climb,V_climb)
Cd_climb = CD(Cl_climb)
Pr = Cd_climb/Cl_climb *W*V_climb

RC = (Hcruise-Hclimbout)/tclimb

P_climb = RC*W + Pr +P_low


#cruise
Cl_cruise = CL(rho_cruise, V_cruise)
Cd_cruise = CD(Cl_cruise)
tcruise = distance/V_cruise

Pr_cruise = Cd_cruise/Cl_cruise *W *V_cruise +P_low

#loiter
Cl_loiter = np.sqrt(3*Cd0*np.pi*A*e)
V_loiter = Velocity(rho_loiter,Cl_loiter)
Cd_loiter = CD(Cl_loiter)
D = Cd_loiter/Cl_loiter * W

P_loiter = D*V_loiter +P_low

#decent
Cl_descent = CL(rho_descent,V_descent)
Cd_descent = CD(Cl_descent)
D = Cd_descent*0.5*rho_descent*V_descent**2*S
T = D + W*np.sin(np.radians(descent_angle))
t_descent = -Hcruise/(V_descent*np.sin(np.radians(descent_angle)))

P_descent = T*V_descent +P_low

#sizing

#battery
power_density = 2000 #W/kg
energy_density = 350 #Wh/kg


inputfactor = 1/(eff_propeller*eff_engine*eff_pmad)

P_fuelcell = Pr_cruise*inputfactor
m_fuelcell = 111.86*np.log(P_fuelcell/1000) - 261.95 + P_fuelcell/3000
Pbat = P_takeoff*inputfactor-P_fuelcell
mbat = Pbat/power_density
Ebat = mbat *energy_density

print("Battery until what phase?: climbout/climb")
batphase = input()

Ebat_needed = (P_takeoff*inputfactor-P_fuelcell)*t_takeoff + (P_climbout*inputfactor-P_fuelcell)*t_climbout
if batphase == "climb":
    Ebat_needed= Ebat_needed + (P_climb*inputfactor-P_fuelcell)*tclimb

if Ebat_needed>Ebat:
    mbat = Ebat_needed/(energy_density*3600)
    Pbat = mbat *power_density
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
E_takeoff = Energy(t_takeoff,P_takeoff,0.45) + Energy(t_climbout,P_climbout,0.45)
E_climb = Energy(tclimb,P_climb,0.50)
E_cruise = Energy(tcruise,Pr_cruise,0.55)
E_descent = Energy(t_descent,P_descent,0.55)
E_loiter = Energy(t_loiter,P_loiter,0.55)
E_taxi = Energy(15*60,Pmax*0.1,0.55)

Etotal = E_takeoff +E_climb+E_cruise+E_descent+E_taxi+E_loiter
Hydrogen = Etotal/E_hydrogen
print(Hydrogen)
print(Etotal/3600000)
print('fight time = ', flight_time/3600)

#radiator
Reynolds = 15000000
T_diff = 50
Pratzl = 0.71
#Nusselt =
print(Pr_cruise,P_loiter)