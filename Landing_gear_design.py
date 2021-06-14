from math import *
from CGloadingdiagram import max_xcg_Hopper_nose

MTOW = 8618
Tip_over = 15                                           #tip over anglee
Scrape_angle = 13                                       #scrape angle
overturn = 55                                           #overturn angle
Phi_clearance = 5

Height_fuselage = 1.91                                  #height fuselage
CGyvalue  = (2/3)*Height_fuselage
height = 0.6
Height_mainlg = CGyvalue+height               #height of cg to ground #estimated at 2/3 of fuselage height + 0.5 landing gear height
wingspan = 20

Tip_height = Height_fuselage                            #estimated with 0.5 meter deflection downwards during taxi operation
Prop_dia = 2.4
Prop_height = Tip_height - Prop_dia/2                   #Estimated to be worst case: with deflection of tip. located at center of the wing
Engine_yloc = 10                                        #lateral engine placement.

n_m = 2                                                 #number of main landing gears
Pn = 12                                               #Percentge of weight on the nose wheel
Lm = sin(radians(Tip_over))*Height_mainlg               #Length of main landing gear to CG
loc = 8.46

dis = loc-(height/tan(radians(13)))


Lm = dis-max_xcg_Hopper_nose
Ln = (100-Pn)/Pn*Lm
print('length of nose landing gear from cg is',Ln)
Pn_total = MTOW*(Pn/100)
Pm_total = MTOW*(1-(Pn/100))/n_m

Lgmain_xloc = Lm+max_xcg_Hopper_nose
Lgnose_xloc = max_xcg_Hopper_nose-Ln
print((Ln**2*tan(radians(overturn))**2)/Height_mainlg-1)
Lgmain_yloc1 = (Ln+Lm)/(sqrt((Ln**2*tan(radians(overturn))**2)/Height_mainlg-1))    #Tip over distance
Lgmain_yloc2 = wingspan/2-Tip_height/tan(radians(Phi_clearance))                    #Tip clearance distance
Lgmain_yloc3 = Engine_yloc-Prop_height/tan(radians(Phi_clearance))                  #Prop clearance distance
#print(Tip_height)
#print(Lgmain_xloc)
print(Lgnose_xloc)                                      #nose wheel location from nose of the aircraft
#print(Lgmain_yloc1,Lgmain_yloc2,Lgmain_yloc3)           #prints the 3 critical y placements from center.



