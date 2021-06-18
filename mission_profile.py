import numpy as np
from ISA import script
import matplotlib.pyplot as plt

W_P = 0.0723
m = 8618 -1814+116
g = 9.81
W = m*g
distance = 2180000# 1600000# 1600000 #480000#1900000 #1880000  # 555600 #900000


S = 43
Cd0 = 0.0251
e = 0.785
A = 9.302

P_low = 14000 #lighting, cockpit, attitude, cooling pump
P_max = W/W_P
print("pmax", P_max)
P_compressor = 64000 #vanaf take-off
P_startupheater = 12000 #taxi
P_charging = 20000
P_charging_decent = 300000

V_climb = 60 #TBD
tclimb = 10*60

#constants
E_hydrogen = 120000000 #J/kg
efficiency_constant = 1.23
eff_engine = 0.93
eff_propeller_climb = 0.8
eff_propeller_cruise = 0.8
eff_pmad = 0.989


rho_climb = 1.05

Hcruise = 3000

#CL_takeoff = 2.4
#rho_takeoff = 1.225
#s_takeoff = 750

rho_cruise = script(Hcruise)
V_cruise = 90

rho_loiter = 1
t_loiter = 45*60 #min

V_descent = 70 #TBD
descent_angle = -1.5 #deg
rho_descent = 1



def CD(Cl):
    return Cd0 + Cl**2/(np.pi*A*e)

def CL(rho,V):
    return W/(0.5*rho*V**2*S)

def Velocity(rho,Cl):
    return np.sqrt(2*W/(rho*S*Cl))

def Energy(t,P, eff_fuelcell):
    return t*P/(eff_fuelcell)

def efficiency(A):
    V = -0.047 * np.log(A) + 0.9782
    return V/efficiency_constant

#startup

P_startup = P_low+P_max*0.05
t_startup = 15*60

#taxi
P_taxi = P_max*0.1 + P_low
t_taxi = 7*60

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

P_engine = P_max/(2*eff_propeller_climb)

P_takeoff = P_max/(eff_propeller_climb*eff_pmad*eff_engine)+P_low +P_compressor
P_climbout = 0.85*P_takeoff
RC_climbout = 4
Hclimbout = t_climbout*RC_climbout

#climb
Cl_climb = CL(rho_climb,V_climb)
Cd_climb = CD(Cl_climb)
Pr = Cd_climb/Cl_climb *W*V_climb
#RC * W
RC = (Hcruise-Hclimbout)/tclimb
#RC = V_climb*(np.sin(np.radians(1.37)))

P_climb = (RC * W + Pr)/(eff_propeller_climb*eff_pmad*eff_engine)+P_compressor+P_low

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
P_engineOEI = P_climbOEI / eff_propeller_climb
print("One engine inoperative", P_climbOEI, P_engineOEI)

#climb one engine inoperative case 1
rho_climbOEI2 = 1.17215
V_climbOEI2 = 1.2 * Velocity(rho_climbOEI2, 1.8)
CL_climbOEI2 = CL(rho_climbOEI2, V_climbOEI2)
Cd_climbOEI2 = CD(CL_climbOEI2)
PrOEI2 = Cd_climbOEI2/CL_climbOEI2 *W*V_climbOEI2
RC_OEI2 = V_climbOEI*(np.sin(np.radians(0.6875)))
P_climbOEI2 = RC_OEI2*W + PrOEI2
P_engineOEI2 = P_climbOEI2 / eff_propeller_climb
print("One engine inoperative case 2", P_climbOEI2, P_engineOEI2)


#cruise
Cl_cruise = CL(rho_cruise, V_cruise)
Cd_cruise = CD(Cl_cruise)
tcruise = distance/V_cruise

Pr_cruise = Cd_cruise/Cl_cruise *W *V_cruise
P_cruise = Pr_cruise/(eff_propeller_cruise*eff_pmad*eff_engine)+P_low+P_compressor+P_charging

#loiter
Cl_loiter = np.sqrt(3*Cd0*np.pi*A*e)
V_loiter = Velocity(rho_loiter,Cl_loiter)
Cd_loiter = CD(Cl_loiter)
D = Cd_loiter/Cl_loiter * W
P_loiter_opt = (D*V_loiter)
P_loiter = (Pr_cruise)/(eff_propeller_cruise*eff_pmad*eff_engine)+P_low+P_compressor

#decent
Cl_descent = CL(rho_descent,V_descent)
Cd_descent = CD(Cl_descent)
D = Cd_descent*0.5*rho_descent*V_descent**2*S
T = D + W*np.sin(np.radians(descent_angle))
t_descent = -Hcruise/(V_descent*np.sin(np.radians(descent_angle)))

P_descent = (T*V_descent)/(eff_propeller_cruise*eff_pmad*eff_engine)+P_low+P_charging_decent+P_compressor
print("descent power", P_descent)
#taxi+shutdown
P_taxishut = P_low+P_max*0.1
t_taxishut = 10*60
#sizing

#battery
power_density = 2000 #W/kg
energy_density = 0.9*369 *3600#Wh/kg



P_fuelcell = P_cruise
m_fuelcell = 111.86*np.log(P_fuelcell/1000) - 261.95 + P_fuelcell/3000
Pbat = P_takeoff-P_fuelcell
mbat = Pbat/power_density
Ebat = mbat * energy_density

print("Battery until what phase?: climbout/climb")
#batphase = input()

Ebat_needed = (P_takeoff-P_fuelcell)*t_takeoff + (P_climbout-P_fuelcell)*t_climbout +P_taxishut*t_taxishut+(P_startup*t_startup)+ P_taxi*t_taxi
Ebat_needed= Ebat_needed + (P_climb-P_fuelcell)*tclimb

if Ebat_needed>Ebat:
    mbat = Ebat_needed/(energy_density)
    Pbat = mbat * power_density
    print("Energy limiting")
else:
    print("power limiting")

print("battery energy used", Ebat_needed/3600000 , "kWh")
print("Fuel cell power = ",P_fuelcell/1000, "kW")
print("Fuel cell mass",m_fuelcell)
print("Battery power = ",Pbat/1000, "kW")
print("Battery energy = ", mbat*energy_density/3600000)
print("Battery mass = ", mbat)

print("battery charged", (P_charging*tcruise+P_charging_decent*t_descent)/3600000)
#Hydrogen calculation
flight_time = t_takeoff+tclimb+tcruise+t_descent+t_loiter
E_takeoff = Energy(t_takeoff,P_takeoff,0.53) + Energy(t_climbout,P_climbout,0.53)
E_climb = Energy(tclimb,P_climb,0.53)
E_cruise = Energy(tcruise,Pr_cruise,0.53)
E_descent = Energy(t_descent,P_descent,0.6)
E_loiter = Energy(t_loiter,P_loiter,0.6)
E_taxi = Energy(t_taxi,P_taxi,0.65)
E_startup = Energy(t_startup,P_startup,0.65)
E_taxishut = Energy(t_taxishut,P_taxishut,0.65)

Etotal = E_takeoff +E_climb+E_cruise+E_descent+E_taxi+E_loiter+E_taxishut+E_startup
Hydrogen = Etotal/E_hydrogen
print("hydrogen", Hydrogen)
print('fight time = ', flight_time/3600)

#radiator
Reynolds = 15000000
T_diff = 50
Pratzl = 0.71
#Nusselt =
print(E_startup/3600000,E_taxi/3600000,E_takeoff/3600000,E_climb/3600000,E_cruise/3600000,E_descent/3600000, E_taxishut/3600000, E_loiter/3600000)
print(P_climb)
print(Etotal/3600000)
print(Cl_cruise/Cd_cruise)
print(RC)