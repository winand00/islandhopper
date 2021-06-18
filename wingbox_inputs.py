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
S_w = 43
taper_wing = 0.5
c_avg_wing = S_w/(b_wing)
c_root_wing= c_avg_wing/((2-taper_wing)/2)
# Root dimensions of the wingbox

w_wing = c_root_wing * 0.5
h_wing = 0.3 * w_wing
# High-lift devices
lx_hld = w_wing / 2 * 1.5  # X-position with respect to the centre of the wingbox
ly_hld = 0.3 * b_wing / 2  # Y-position with respect to the centre of the wingbox
delta_CL = 1.1
#F_hld = 15000     # Force of hld
F_hld = 0.5*delta_CL*1.225*S_w*V_F**2 / 2

# Aileron
lx_a = w_wing / 2  # X-position with respect to the centre of the wingbox
ly_a = 0.8 * b_wing / 2  # Y-position with respect to the centre of the wingbox
F_a = 0.5 * 1.225 * V_cruise ** 2 * S_w * 0.25 / w_wing  # Force of aileron

# Vertical tail
#S_vert = 10 * 1.1
S_vert = 9.7
A_vert = 2.2
b_vert = np.sqrt(A_vert*S_vert) * 2

taper_vert = 0.6
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
F_r = 5115

# Horizontal tail

taper_hor = 0.8
S_h = 8
A_h = 5.193
b_hor = np.sqrt(S_h*A_h)
c_avg_hor = S_h/(b_hor)
c_root_hor = c_avg_hor/((1-taper_hor)/2+taper_hor)
print(c_root_hor)
CL_h = 0.35*A_h**(1/3)
L_hor = -0.5*CL_h*1.225*S_h*V_cruise**2 / b_hor  # Lift distribution on the horizontal tail with n = 1

# Root dimensions of the wingbox
w_hor = 0.5 * c_root_hor
h_hor = 0.3 * w_hor

# Elevator
ly_el = b_hor / 2 / 2  # y position of the elevator
lx_el = w_hor / 2 * 1.5  # x position of the elevator
F_el = tail_load_elevator(lf_pos)  # elevator force


#Systems inputs:
#Engine
w_engine = 161 * 9.81
ly_e = 0.75 * b_wing / 2

lz_e = h_wing / 2
powersetting = 0.5
w_radiator = 95 * 9.81
w_prop = 79 * 9.81

#Fuel cell
w_fc = (117 + 50 + 137.3) * 9.81
ly_fc = 0.33 * b_wing

#Battery
w_bat = 210 * 9.81
ly_bat = 0.66 * b_wing

#Other systems
w_sys = w_radiator
ly_sys = 0.7 * b_wing / 2

# Structures
l_fuselage = 12.2 + 1
x_pos_wing = 5.22
D_fuselage = 2.25
weight_w = (412.02 * 2 + 150) * 9.81
w_systems = w_engine+w_radiator+w_prop+w_bat+w_sys
weight_wing = (w_systems+weight_w)