from math import *
from CGloadingdiagram import max_xcg_Hopper_nose
print()
print()

##inputs
height = 0.1                                            #Landing gear height from bottom fuselage.
Height_fuselage = 1.91                                  #height fuselage
CGyvalue= (2/3)*Height_fuselage                         #y valeu of cg from bottom of the fuselage. (estimated at 2/3)
Engine_yloc = 5                                         #lateral engine placement.
n_m = 2                                                 #number of main landing gears
Pn = 5                                                  #Percentge of weight on the nose wheel
minnoselg_xpos = 0.5                                    #minimum distance for the nose lg to nose of aircraft

MTOW = 8618
wingspan = 20
Prop_dia = 2.4
loc = 8.46                                              #Beginning of tailcone from nose

Tip_over = 15                                           #tip over anglee
Scrape_angle = 13                                       #scrape angle
overturn = 55                                           #overturn angle
Phi_clearance = 5

## calculations
Tip_height = Height_fuselage+height-0.5                 #estimated with 0.5 meter deflection downwards during taxi operation
Prop_height = Tip_height - Prop_dia/2                   #Estimated to be worst case: with deflection of tip. located at center of the wing
Height_CG = CGyvalue+height                             #height of cg to ground
Lm1 = Height_CG*tan(radians(Tip_over))                  #Length of main landing gear to CG (for tip_over
dis = loc-(height/tan(radians(13)))                     #x location of main lg from nose (for scrape)
Lm2 = dis-max_xcg_Hopper_nose

##Longititudinal placement iteration for scrape and tip over angle

while abs(Lm1-Lm2)> 0.001:
    height = height+0.0001
    Height_CG = CGyvalue + height                       # height of cg to ground
    Lm1 = Height_CG * tan(radians(Tip_over))            # Length of main landing gear to CG (for tip_over
    dis = loc - (height / tan(radians(13)))             # x location of main lg from nose (for scrape)
    Lm2 = dis - max_xcg_Hopper_nose

Lm=Lm1
Ln = (100-Pn)/Pn*Lm
Pn_total = MTOW*(Pn/100)
Pm_total = MTOW*(1-(Pn/100))/n_m
Lgmain_xloc = Lm+max_xcg_Hopper_nose
Lgnose_xloc = max_xcg_Hopper_nose-Ln

while(Lgnose_xloc<minnoselg_xpos or Pn<8):
    Pn=Pn+0.001
    Ln = (100 - Pn) / Pn * Lm
    Lgnose_xloc = max_xcg_Hopper_nose - Ln


print('Nose landing gear from distance from nose is',round(Lgnose_xloc,3),'[m], and the load per landing gear is:',round(Pn_total,3),'[N]')
print('Main landing gear from distance from nose is',round(Lgmain_xloc,3),'[m], and the load per landing gear is:',round(Pm_total,3),'[N]')

##Lateral positioning


Lgmain_yloc1 = (Ln+Lm)/(sqrt((Ln**2*tan(radians(overturn))**2)/Height_CG-1))    #Tip over distance
Lgmain_yloc2 = wingspan/2-Tip_height/tan(radians(Phi_clearance))                    #Tip clearance distance
Lgmain_yloc3 = Engine_yloc-Prop_height/tan(radians(Phi_clearance))                  #Prop clearance distance

if max(Lgmain_yloc1,Lgmain_yloc2,Lgmain_yloc3) == Lgmain_yloc1:
    print('Tip over is critical. Distance from center is:               ',round(Lgmain_yloc1,2),'[m]')
elif max(Lgmain_yloc1,Lgmain_yloc2,Lgmain_yloc3) == Lgmain_yloc2:
    print('Tip clearance is critical. Distance from center is:          ',round(Lgmain_yloc2,2),'[m]')
elif max(Lgmain_yloc1,Lgmain_yloc2,Lgmain_yloc3) == Lgmain_yloc3:
    print('Prop clearance is critical. Distance from center is:         ',round(Lgmain_yloc3,2),'[m]')
print('Load percentage of nose wheel is:                            ',round(Pn,3),'[%]')
print('Height of the landing gear is:                               ',round(height,3),'[m]')
#print(Lgmain_yloc1,Lgmain_yloc2,Lgmain_yloc3)          #prints the 3 critical y placements from center.




