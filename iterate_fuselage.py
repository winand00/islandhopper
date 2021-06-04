import numpy as np
from fuselage import create_fuselage


def iterate_fuselage(t_skin, n_str):
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for str_n in n_str:
            fuselage = create_fuselage(t, str_n)
            if not fuselage.is_failing():
                working_designs[f't_sk:{round(t, 3)}, n_str:{str_n}'] = round(fuselage.weight, 2)
            else:
                failed_designs[f't_sk:{round(t, 3)}, n_str:{str_n}'] = fuselage.is_failing()
    min_weight = min(working_designs.values())
    return min_weight, working_designs, failed_designs



n_str = np.arange(12, 40, 4)
t_skin = np.arange(0.002, 0.006, 0.001)

print(iterate_fuselage(t_skin, n_str))