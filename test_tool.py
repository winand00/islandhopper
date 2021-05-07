import unittest
import tool

class TestTool(unittest.TestCase):

    def test_calculate_cl_opt(self):
        result = tool.calculate_cl_opt(0.015, 10, 0.7)
        self.assertEqual(round(result, 5), 0.18708)
        self.assertTrue(result > 0)

    def test_calculate_v_cruise(self):
        result = tool.calculate_v_cruise(6950 * 9.81, 25.08, 0.54895, tool.calculate_cl_opt(0.015, 10, 0.7))
        self.assertEqual(round(result, 5), 230.08840)
        self.assertTrue(result > 0)

if __name__ == "__main__":
    unittest.main()

