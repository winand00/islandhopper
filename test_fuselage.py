from fuselage import Fuselage, Stringer
import unittest

t_sk = 0.005
n_str = 4


def create_fuselage(t_sk, n_str):
    weight_ac = 84516
    n = 1
    length = 10
    weight_wing = 10000
    weight_dist = (weight_ac - weight_wing) / length * n
    t = t_sk
    D = 1
    x_t = length - 0.1
    F_t = 1000
    x_w = 5
    F_w = 0  # weight_ac * n - weight_wing + F_t
    T_t = 500
    buckling_factor = 0  # 5 # fraction of sigma_y
    # Material
    # Stringer
    rho_str = 2820
    E_str = 69 * 10 ** 9
    sigma_y_str = 450 * 10 ** 6
    poisson_str = 0.33

    # Skin
    rho_sk = 2820
    E_sk = 69 * 10 ** 9
    sigma_y_sk = 450 * 10 ** 6
    poisson_sk = 0.33

    # G = 26.4 * 10 ** 9

    # Make stringer
    t_str = 0.005
    w_str = 0.1
    h_str = 0.1

    stringer = Stringer(t_str, w_str, h_str, rho_str, sigma_y_str, E_str, poisson_str)
    n_str = n_str  # Has to be a multiple of 4

    fuselage = Fuselage(stringer, n_str, length, t, D, weight_dist, x_w, T_t, x_t, F_t, rho_sk, E_sk, sigma_y_sk,
                        poisson_sk, buckling_factor)
    return fuselage

fuselage = create_fuselage(t_sk, n_str)

class StressTest(unittest.TestCase):
    def test_skin_area(self):




if __name__ == "__main__":
        unittest.main()