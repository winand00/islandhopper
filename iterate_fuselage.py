import numpy as np
from fuselage import create_fuselage
from create_wingbox import AL, AL7040


def iterate_fuselage(t_skin, n_str, material_skin, material_stringer):
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for str_n in n_str:
            fuselage = create_fuselage(t, str_n, material_skin, material_stringer)
            if not fuselage.is_failing():
                working_designs[f't_sk:{round(t, 3)}, n_str:{str_n}'] = round(fuselage.weight, 2)
            else:
                failed_designs[f't_sk:{round(t, 3)}, n_str:{str_n}'] = fuselage.is_failing()
    min_weight = min(working_designs.values())
    return min_weight, working_designs, failed_designs


material_skin = AL
material_stringer = AL7040
n_str = np.arange(12, 40, 4)
t_skin = np.arange(0.002, 0.006, 0.001)

print(iterate_fuselage(t_skin, n_str, material_skin, material_stringer))