import unittest
from hydrogen_tank_v2 import Tank

import numpy as np
import math

# CONSTANTS
k_overlap = 0.75      #Used in thesis
P = 5.7
# sigma_allowable = 1

t_ply = 0.00014       #[m] #Table 4-3
t_metal_liner = 0.001
t_polyamide_liner = 0.001

rho_comp = 2700
k_comp = 60
rho_ins = 200
rho_metal_liner = 2700
rho_polyamide_liner = 1010
k_ins = 0.025

y = 1.4
R_over_t = 130
s_over_R = 0.35           #(fig 4.22)
r = 0.1                   #(fig 4.23a) [R_fillet/R]
t_junction_over_R = 0.027 #(fig 4.23d)
phi = 24                  #degrees (fig 4.23b) [layup angle]
i_ratio = 2.5             #(fig 4.23b)   [ratio between 0 degrees fibers and phi degrees fibres]
k_overlap = 0.75          #(fig 5.5)
BOR = 0.016

rho_hydr = 70 # [kg/m3]

d_H_vap = 446100 #J/kg
BOR_percentage = 0.016 # [-]
d_T = 293.15 # temp diff between 20 K and 40 C
#h_out = -0.23083#30 #[W/m2K]

class TestHydrogenTank(unittest.TestCase):

    def setUp(self):
        # self.tank_1 = Tank(20, 20, 20, 0.2, 1) # Many spheres
        # self.tank_2 = Tank(2, 2, 2, 0.01, 1)   # Small radius
        # self.tank_3 = Tank(2, 2, 2, 10, 1)      # Large radius
        # self.tank_4 = Tank(1, 1, 1, 1, 1)      # Sphere
        # self.tank_5 = Tank(1, 1, 10, 1, 1)      # Long line


    # def test_multi_cell_dimensions(self):
    #     # m, n, p, R, t_ins, t_polyamide_liner
    #     # t_polyamide_liner = 0.001
    #     result = hydrogen_tank_v2.multi_cell_dimensions(2, 2, 2, 100, 1, 1)
    #     self.assertTrue(result[0] == result[1] == result[2])
    #
    # def test_multi_cell_v(self):
    #     # (m, n, p, R)
    #     # result = hydrogen_tank_v2.multi_cell_v()
    #     pass
    #
    # def test_multi_cell_s(self):
    #     # (m, n, p, R)
    #     # result = hydrogen_tank_v2.multi_cell_s()
    #     pass
    #
    # def test_multi_cell_m_incl_insulation(self):
    #     # (m, n, p, R, Mass)
    #     # result = hydrogen_tank_v2.multi_cell_m_incl_insulation()
    #     pass


if __name__ == "__main__":
    unittest.main()