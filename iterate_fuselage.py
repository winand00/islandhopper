import numpy as np
from fuselage import create_fuselage
from create_wingbox import AL6061, Glare, AL2099, AL7055, AL7040


def iterate_fuselage(t_skin, n_str, material_skin, materials_stringer):
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for str_n in n_str:
            for material_stringer in materials_stringer:
                fuselage = create_fuselage(t, str_n, material_skin, material_stringer)
                if not fuselage.is_failing():
                    working_designs[f't_sk:{round(t, 4)}, n_str:{str_n}, Material stringer:{material_stringer.name}'] = round(fuselage.weight, 2)
                else:
                    failed_designs[f't_sk:{round(t, 4)}, n_str:{str_n}, Material stringer:{material_stringer.name}'] = fuselage.is_failing()
    min_weight = min(working_designs.values())
    return min_weight, min(working_designs)#, failed_designs


material_skin = AL6061
material_stringer = [AL7040, AL2099, AL7055]
n_str = np.arange(4, 40, 4)
t_skin = np.arange(0.002, 0.006, 0.0005)

print(iterate_fuselage(t_skin, n_str, material_skin, material_stringer))