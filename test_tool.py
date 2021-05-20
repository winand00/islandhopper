import unittest
import tool

# Inputs jetstream 31

cd_0 = 0.0376               # Zero lift drag
A = 10                      # Aspect ratio
e = 0.75                    # Oswald efficiency
W = 68179.5                 # Weight of aircraft [N]
rho = 0.54895               # Density at cruise h [kg/m3]
rho_sealevel = 1.225        # Density at sealevel [kg/m3]
S = 25.08                   # Wing  [m2]
specific_energy = 46200000  # Specific energy of fuel [J/kg]
m_energy = 1483             # Mass of fuel [kg]
m = 6950                    # Mass of aircraft [kg]
L_over_D = 8.16             # Lift over drag ratio
efficiency_fuelcell = 0.9   # Efficiency fuel cell
efficiency_prop = 0.8       # Efficiency propeller
p_max = 1342500             # Max power [W]
cl_takeoff = 1.8            # Lift coefficient at take-off
fuel_W_maxpayload = 664.1   # Fuel weight @ max payload
cl_max = 2                  # Maximum lift coefficient
D = 2.69                    # Propeller diameter
B = 4                       # Number of blades per propeller
N = 2                       # Number of propellers
efficiency_r = 0.1          # Fraction of total used energy that is recovered for other systems
fuel_weight_maxfuel = 1483  # Fuel weight
fuel_speed = 1000 / 60      # Fuel speed in l/min
battery = False             # Aircraft on batteries

class TestTool(unittest.TestCase):

    def test_calculate_cl_opt(self):
        result = tool.calculate_cl_opt(cd_0, A, e)
        self.assertTrue(result > 0)
        self.assertEqual(round(result, 5), 0.94124)

    def test_calculate_v_cruise(self):
        result = tool.calculate_v_cruise(W, S, rho, tool.calculate_cl_opt(cd_0, A, e))
        self.assertTrue(result > 0)
        self.assertEqual(round(result, 5), 102.57984)

    def test_calculate_range(self):
        result = tool.calculate_range(specific_energy, fuel_W_maxpayload, m, L_over_D, 1, efficiency_prop)
        self.assertTrue(result > 0)
        self.assertEqual(round(result, 3), 2937661.757)

    def test_calculate_max_climb_rate_and_gradient(self):
        result = tool.calculate_max_climb_rate_and_gradient(p_max, W, S, cd_0, rho, A, e, efficiency_prop)
        self.assertEqual(round(result[0], 5), 8.56186)
        self.assertEqual(round(result[1], 5), 0.10985)

    def test_calculate_takeoff_parameter(self):
        result = tool.calculate_takeoff_parameter(W, S, p_max, cl_takeoff)
        self.assertEqual(round(result, 5), 76.69963)

    def test_calculate_landing_distance(self):
        result = tool.calculate_landing_distance(W, rho_sealevel, S, cl_max)
        self.assertTrue(0 < result < 5500)
        self.assertEqual(round(result, 5), 1312.6379)

    def test_calculate_refuel_time(self):
        result = tool.calculate_refuel_time(fuel_W_maxpayload, fuel_speed)
        self.assertTrue(result > 0)
        self.assertEqual(round(result, 5), 39.846)

    def test_calculate_charging_time(self):
        result = tool.calculate_charging_time(650, 250, 0.95)
        self.assertTrue(result > 0)
        self.assertEqual(round(result, 2), 2.74)

    def test_calculate_max_sound_pressure_level(self):
        result = tool.calculate_max_sound_pressure_level(D, B, N, p_max, efficiency_prop)
        self.assertTrue(0 < result < 150)
        self.assertEqual(round(result, 2), 81.76)

if __name__ == "__main__":
    unittest.main()

