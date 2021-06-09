from shear_stress import Stress

import unittest


h = 1
w = 1
Izz = 1
Ixx = 1
t = 0.005
Vx = 1000
Vz = 1000
Mx = 500
Mz = 500
T = 50
x_centroid = 0.5
z_centroid = 0.5
sigma_F = 0

stress = Stress(h, w, Izz, Ixx, t, Vx, Vz, Mx, Mz, T, x_centroid, z_centroid, sigma_F)

class StressTest(unittest.TestCase):
    def test_bending_stress_x(self):
        sigma = stress.bendingstress_x(1,1)
        self.assertAlmostEqual(sigma, 250, 3)
        sigma = stress.bendingstress_x(1,0.5)
        self.assertAlmostEqual(sigma, 0, 3)

    def test_bending_stress_z(self):
        sigma = stress.bendingstress_z(1,1)
        self.assertAlmostEqual(sigma, -250, 3)
        sigma = stress.bendingstress_z(1/2,1)
        self.assertAlmostEqual(sigma, 0, 3)

    def test_bending_stress_xz(self):
        sigma = stress.bendingstress_xz(1,1)
        self.assertAlmostEqual(sigma, 500, 3)
        sigma = stress.bendingstress_xz(0,1)
        self.assertAlmostEqual(sigma, 0, 3)

    def test_shear_stress_x(self):
        tau = stress.shearstress_x(0,0.5)
        self.assertAlmostEqual(tau, 0, 3)
        tau = stress.shearstress_x(0, 1)
        self.assertAlmostEqual(tau, -250, 3)
        tau = stress.shearstress_x(0, 0)
        self.assertAlmostEqual(tau, 250, 3)
        tau = stress.shearstress_x(0.5, 1)
        self.assertAlmostEqual(tau, -375, 3)

    def test_shear_stress_z(self):
        tau = stress.shearstress_z(0.5, 1)
        self.assertAlmostEqual(tau, 0, 3)
        tau = stress.shearstress_z(0, 0)
        self.assertAlmostEqual(tau, -250, 3)
        tau = stress.shearstress_z(0, 0.5)
        self.assertAlmostEqual(tau, -375, 3)

    def test_shear_torsion(self):
        tau = stress.shear_torque(0.5, 1)
        self.assertAlmostEqual(tau, 5000, 3)
        tau = stress.shear_torque(1, 0.5)
        self.assertAlmostEqual(tau, 5000, 3)

    def test_shear_total(self):
        tau = stress.shear_total(0, 0)
        self.assertAlmostEqual(tau, 5000, 3)
        tau = stress.shear_total(0, 0.5)
        self.assertAlmostEqual(tau, 4625, 3)
        tau = stress.shear_total(0.5, 1)
        self.assertAlmostEqual(tau, 4625, 3)
        tau = stress.shear_total(1, 1)
        self.assertAlmostEqual(tau, 5000, 3)

    def test_von_mises(self):
        sigma = stress.von_mises(1, 1)
        self.assertAlmostEqual(sigma, 8674.676, 3)



if __name__ == "__main__":


    unittest.main()