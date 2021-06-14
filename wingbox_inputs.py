# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:56:14 2021

@author: DickD_000
"""
import numpy as np
from design_loads import tail_load_elevator, design_loads
lf_pos = design_loads()[0]
# Aerodynamic inputs:
L_D = 9.64
V_cruise = 90
V_F = 69.11

# Main wing
b_wing = 20
S_w = 40
taper_wing = 0.5

# Root dimensions of the wingbox
h_wing = 0.35
w_wing = 1.33 #1.5
# High-lift devices
lx_hld = w_wing / 2 * 1.5  # X-position with respect to the centre of the wingbox
ly_hld = 0.3 * b_wing / 2  # Y-position with respect to the centre of the wingbox
delta_CL = 0.78
#F_hld = 15000     # Force of hld
F_hld = 0.5*delta_CL*1.225*S_w*V_F**2 / 2

# Aileron
lx_a = w_wing / 2  # X-position with respect to the centre of the wingbox
ly_a = 0.8 * b_wing / 2  # Y-position with respect to the centre of the wingbox
F_a = 0.5 * 1.225 * V_cruise ** 2 * S_w * 0.25 / w_wing  # Force of aileron

# Vertical tail
#S_vert = 10 * 1.1
S_vert = 16.1 * 1.1
A_vert = 1.2
b_vert = np.sqrt(A_vert*S_vert) * 2
taper_vert = 0.8
c_avg_vert = S_vert/(b_vert/2)
c_root_vert = c_avg_vert/((1-taper_vert)/2+taper_vert)
#L_ver = 5300 / (b_vert/2)
L_ver = 14600 / (b_vert/2)  # Lift distribution on the vertical tail
lz_hor = b_vert/4  # Z-position of the horizontal tail on the vertical tail
# Root dimensions of the wingbox
w_vert = 0.25 * c_root_vert
h_vert = w_vert * 0.3

#Rudder
ly_r = b_vert / 2 / 2
lx_r = 3 * w_vert / 2
#F_r = 3000
F_r = 9300

# Horizontal tail

taper_hor = 0.8
S_h = 11.6
A_h = 6.93
b_hor = np.sqrt(S_h*A_h)
c_avg_hor = S_h/(b_hor)
c_root_hor = c_avg_hor/((1-taper_hor)/2+taper_hor)
CL_h = 0.35*A_h**(1/3)
L_hor = -0.5*CL_h*1.225*S_h*V_cruise**2 / b_hor  # Lift distribution on the horizontal tail with n = 1
print(L_hor)
# Root dimensions of the wingbox
w_hor = 0.5 * c_root_hor
h_hor = 0.3 * w_hor

# Elevator
ly_el = b_hor / 2 / 2  # y position of the elevator
lx_el = w_hor / 2 * 1.5  # x position of the elevator
F_el = tail_load_elevator(lf_pos)  # elevator force


#Systems inputs:
#Engine
w_engine = 180 * 9.81
ly_e = 1 * b_wing / 2

lz_e = h_wing / 2
powersetting = 0.5
w_radiator = 133 * 9.81
w_prop = 81 * 9.81

#Fuel cell
w_fc = (117 + 50 + 142.8) * 9.81
ly_fc = 0.5 * b_wing

#Battery
w_bat = 222 * 9.81
ly_bat = 0.4 * b_wing

#Other systems
w_sys = 113 * 9.81 + w_radiator
ly_sys = 0.3 * b_wing

# Structures
l_fuselage = 12.2
x_pos_wing = l_fuselage / 2
D_fuselage = 2.09
weight_wing = 1000 * 9.81