from create_wingbox import make_wingbox, AL7040, AL, Ti, Glare, n_ult_pos, n_ult_neg
import numpy as np


def iterate_wing_box(t_skin, n_str, str_size, material, n, type):
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for str_n in n_str:
            for size in str_size:
                wingbox_pos = make_wingbox(t, str_n, size, material, n[0], type)
                wingbox_neg = make_wingbox(t, str_n, size, material, n[1], type)
                if not wingbox_pos.is_failing()[0] and not wingbox_neg.is_failing()[0]:
                    working_designs[f't_sk:{round(t, 4)}, n_str:{str_n*5}, size_str:{size}'] = round(wingbox_pos.weight, 2)
                else:
                    failed_designs[f't_sk:{round(t, 4)}, n_str:{str_n*5}, size_str:{size}'] = wingbox_pos.is_failing()[1]
    min_weight = min(working_designs.values())
    return min_weight, min(working_designs), failed_designs

type = 'wing'
size_str = np.arange(0.03, 0.04, 0.01)
n_str = np.arange(0, 8, 1)
t_skin = np.arange(0.0035, 0.005, 0.0005)
n = [n_ult_pos, n_ult_neg]
print(iterate_wing_box(t_skin, n_str, size_str, AL7040, n, type))
#print(iterate_wing_box(t_skin, n_str, size_str, AL, n, type))
#print(iterate_wing_box(t_skin, n_str, size_str, Ti, n, type))
#print(iterate_wing_box(t_skin, n_str, size_str, Glare, n, type))