from create_wingbox import make_wingbox, AL7040, AL7055, AL2099, AL, Ti, Glare, n_ult_pos, n_ult_neg, n_ult_flaps, n_ult_ail
import numpy as np


def iterate_wing_box(t_skin, n_str_top, n_str_bot, str_size, materials, n, type):
    lowest_weight = np.float('inf')
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for top_str_n in n_str_top[::-1]:
            for bot_str_n in n_str_bot[::-1]:
                for size in str_size[::-1]:
                    for skin_material in materials:
                        for top_str_material in materials:
                            for bot_str_material in materials:
                                wingbox_pos = make_wingbox(t, top_str_n, bot_str_n, size, skin_material, top_str_material, bot_str_material, n[0], type)
                                pos_fail = wingbox_pos.is_failing()[0]
                                if pos_fail or wingbox_pos.weight > lowest_weight:
                                    break
                                wingbox_neg = make_wingbox(t, top_str_n, bot_str_n, size,  skin_material, top_str_material, bot_str_material, n[1], type)
                                neg_fail = wingbox_neg.is_failing()[0]
                                if neg_fail:
                                    break
                                wingbox_flaps = make_wingbox(t, top_str_n, bot_str_n, size,  skin_material, top_str_material, bot_str_material, n[2], type)
                                flap_fail = wingbox_flaps.is_failing()[0]
                                if flap_fail:
                                    break
                                wingbox_ail = make_wingbox(t, top_str_n, bot_str_n, size,  skin_material, top_str_material, bot_str_material, n[3], type)
                                ail_fail = wingbox_ail.is_failing()[0]
                                if pos_fail:
                                    break
                                if not pos_fail and not neg_fail and not flap_fail and not ail_fail:
                                    lowest_weight = wingbox_pos.weight
                                    working_designs[f't_sk:{round(t, 4)}, n_str_top:{top_str_n*5}, n_str_bot:{bot_str_n*2}, size_str:{size}, skin:{skin_material.name}, top_str:{top_str_material.name}, bot_str:{bot_str_material.name}'] = round(wingbox_pos.weight, 2)
                                else:
                                    failed_designs[f't_sk:{round(t, 4)}, n_str_top:{top_str_n*5}, n_str_bot:{bot_str_n*2}, size_str:{size}, skin:{skin_material.name}, top_str:{top_str_material.name}, bot_str:{bot_str_material.name}'] = wingbox_pos.is_failing()[1]

    min_weight = min(working_designs.values())
    return min_weight, min(working_designs)#, failed_designs

type = 'wing'
size_str = np.arange(0.03, 0.05, 0.01)
n_str_top = np.arange(2, 6, 1)
n_str_bot = np.arange(1, 5, 1)
t_skin = np.arange(0.001, 0.003, 0.0005)
n = [n_ult_pos, n_ult_neg, n_ult_flaps, n_ult_ail]
materials = [AL7055, AL2099]
print(iterate_wing_box(t_skin, n_str_top, n_str_bot, size_str, materials, n, type))
#print(iterate_wing_box(t_skin, n_str, size_str, AL, n, type))
#print(iterate_wing_box(t_skin, n_str, size_str, Ti, n, type))
#print(iterate_wing_box(t_skin, n_str, size_str, Glare, n, type))