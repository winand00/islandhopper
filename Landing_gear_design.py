from math import *
from CGloadingdiagram import max_xcg_Hopper_nose

MTOW = 8618
Tip_over = 15           #tip over angle
Scrape_angle = 13       #scrape angle
overturn = 55           #overturn angle
Phi_clearance = 5

Height_fuselage = 1.91  #height fuselage
Height_mainlg = (2/3)*Height_fuselage+0.5
wingspan = 20

Tip_height = Height_fuselage
Prop_dia = 2.4
Prop_height = Tip_height - Prop_dia/2
Engine_yloc = 10
n_m = 2                 #number of main landing gears
Pn = 8                  #Percentge of weight on the nose wheel
Lm = sin(radians(Tip_over))*Height_mainlg                  #Length of main landing gear to CG


Ln = (100-Pn)/Pn*Lm
Pn_total = MTOW*(Pn/100)
Pm_total = MTOW*(1-(Pn/100))/n_m

Lgmain_xloc = Lm+max_xcg_Hopper_nose
Lgnose_xloc = max_xcg_Hopper_nose-Ln
Lgmain_yloc1 = (Ln+Lm)/(sqrt((Ln**2*tan(radians(overturn))**2)/Height_mainlg-1))    #Tip over distance
Lgmain_yloc2 = wingspan/2-Tip_height/tan(radians(Phi_clearance))                    #Tip clearance distance
Lgmain_yloc3 = Engine_yloc-Prop_height/tan(radians(Phi_clearance))                  #Prop clearance distance
print(Tip_height)
print(Lgmain_xloc)
print(Lgnose_xloc)
print(Lgmain_yloc1,Lgmain_yloc2,Lgmain_yloc3)



