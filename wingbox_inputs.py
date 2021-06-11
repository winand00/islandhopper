# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:56:14 2021

@author: DickD_000
"""
import numpy as np
# Aerodynamic inputs:
L_D = 9.64
V_cruise = 90


# Main wing
b_wing = 20
S_w = 45
taper_wing = 0.5

# Root dimensions of the wingbox
h_wing = 0.35
w_wing = 1.33 #1.5
# High-lift devices
lx_hld = w_wing / 2  # X-position with respect to the centre of the wingbox
ly_hld = 0.3 * b_wing / 2  # Y-position with respect to the centre of the wingbox
F_hld = 15000     # Force of hld
# Aileron
lx_a = w_wing / 2  # X-position with respect to the centre of the wingbox
ly_a = 0.8 * b_wing / 2  # Y-position with respect to the centre of the wingbox
F_a = 0.5 * 1.225 * V_cruise ** 2 * S_w * 0.25 / w_wing  # Force of aileron

# Vertical tail
S_vert = 15
A_vert = 1
b_vert = np.sqrt(A_vert*S_vert) * 2
taper_vert = 0.8
c_avg_vert = S_vert/b_vert
c_root_vert = c_avg_vert/((1-taper_vert)/2+taper_vert)
L_ver = 17500 / b_vert  # Lift distribution on the vertical tail
lz_hor = b_vert/3  # Z-position of the horizontal tail on the vertical tail
# Root dimensions of the wingbox
h_vert = 0.4 * c_root_vert
w_vert = h_vert * 0.3
#Rudder
ly_r = b_vert / 2 / 2
lx_r = w_vert / 2
F_r = 1000

# Horizontal tail
b_hor = 6.736 # From L410
taper_hor = 0.8
L_hor = 3500 / b_hor  # Lift distribution on the horizontal tail with n = 1
# Root dimensions of the wingbox
h_hor = 0.1
w_hor = 0.5
# Elevator
ly_el = b_hor / 2 / 2  # y position of the elevator
lx_el = w_hor / 2  # x position of the elevator
F_el = 1000  # elevator force


#Systems inputs:
#Engine
w_engine = 130
ly_e = 1 * b_wing / 2
lz_e = h_wing / 2
powersetting = 0.5
w_radiator = 100
w_prop = 80

#Fuel cell
w_fc = 140
ly_fc = ly_e

#Battery
w_bat = 135
ly_bat = ly_e

# Structures
l_fuselage = 10
x_pos_wing = l_fuselage / 2
D_fuselage = 1.91
weight_wing = 1000 * 9.81