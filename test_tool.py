import unittest
import tool

class TestTool(unittest.TestCase):

    def test_calculate_cl_opt(self):
        result = tool.calculate_cl_opt(0.015, 10, 0.7)
        self.assertEqual(round(result, 5), 0.18708)

    def test_calculate_v_cruise(self):
        result = tool.calculate_v_cruise(7350 * 9.81, 25.08, 0.54895, tool.calculate_cl_opt(0.015, 10, 0.7))
        self.assertEqual(round(result, 5), 236.61703)

    def test_calculate_range(self):
        result = tool.calculate_range(46200000, 664.1, 7350, 8.16, 0.65)
        self.assertEqual(round(result, 3), 2256953.569)

    def test_calculate_max_climb_rate_and_gradient(self):
        result = tool.calculate_max_climb_rate_and_gradient(1074000, 7350 * 9.81, 25.08, 0.0376, 0.54895, 10.0, 0.75)
        self.assertEqual(round(result[0], 5), 7.50055)
        self.assertEqual(round(result[1], 5), 0.09358)

    def test_calculate_takeoff_parameter(self):
        result = tool.calculate_takeoff_parameter(7350 * 9.81, 25.08, 1074000, 1.8)
        self.assertEqual(round(result, 5), 107.22804)

    def test_calculate_runway_length_landing(self):
        result = tool.calculate_runway_length_landing(7350 * 9.81, 0.54895, 10.0, 2)
        self.assertEqual(round(result, 5), 3107694.34375)
        self.assertTrue(0 < result < 5500000)

    def test_calculate_refuel_time(self):
        result = tool.calculate_refuel_time(1483, 1000 / 60)
        self.assertEqual(round(result, 5), 88.98)

    def test_calculate_charging_time(self):
        result = tool.calculate_charging_time(650, 250, 0.95)
        self.assertEqual(round(result, 2), 2.74)

    def test_calculate_max_sound_pressure_level(self):
        result = tool.calculate_max_sound_pressure_level(1342500, 2.49936, 0.9, 5, 2, 50)
        self.assertEqual(round(result, 2), 163.88)

    def test_



if __name__ == "__main__":
    unittest.main()

