from wing_box import make_wingbox
import numpy as np


def iterate_wing_box(t_skin, n_str, size_str):
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for n in n_str:
            for size in size_str:
                wingbox = make_wingbox(t, n, size)
                if not wingbox.is_failing()[0]:
                    working_designs[f't_sk:{round(t, 3)}, n_str:{n*5}, size_str:{size}'] = round(wingbox.weight, 2)
                else:
                    failed_designs[f't_sk:{round(t, 3)}, n_str:{n*5}, size_str:{size}'] = wingbox.is_failing()[1]
    return working_designs, failed_designs


size_str = np.arange(0.03, 0.06, 0.01)
n_str = np.arange(0, 5, 1)
t_skin = np.arange(0.003, 0.005, 0.001)
print(iterate_wing_box(t_skin, n_str, size_str))
print(min(iterate_wing_box(t_skin, n_str, size_str)[0].values()))