from create_wingbox import make_wingbox, AL7040, AL7055, AL2099, AL, Ti, Glare, n_ult_pos, n_ult_neg, n_ult_flaps, n_ult_ail
import numpy as np
from time import time


def iterate_wing_box(t_skin, n_str_top, n_str_bot, str_size, materials, n, type):
    lowest_weight = float('inf')
    working_designs = {}
    failed_designs = {}
    for t in t_skin:
        for top_str_n in n_str_top:
            for bot_str_n in n_str_bot:
                for size in str_size:
                    for skin_material in materials[1:]:
                        for top_str_material in materials[:1]:
                            for bot_str_material in materials[1:]:
                                wingbox_pos = make_wingbox(t, top_str_n, bot_str_n, size, skin_material, top_str_material, bot_str_material, n[0], type)
                                pos_fail = wingbox_pos.is_failing()[0]
                                wb_weight = wingbox_pos.weight
                                if pos_fail or wb_weight > lowest_weight:
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
                                    lowest_weight = wb_weight
                                    working_designs[f't_sk:{round(t, 5)}, n_str_top:{top_str_n*5}, n_str_bot:{bot_str_n*2}, size_str:{size}, skin:{skin_material.name}, top_str:{top_str_material.name}, bot_str:{bot_str_material.name}'] = round(wb_weight, 2)
                                else:
                                    failed_designs[f't_sk:{round(t, 5)}, n_str_top:{top_str_n*5}, n_str_bot:{bot_str_n*2}, size_str:{size}, skin:{skin_material.name}, top_str:{top_str_material.name}, bot_str:{bot_str_material.name}'] = wingbox_pos.is_failing()[1]

    min_weight = min(working_designs.values())
    return min_weight, lowest_weight, min(working_designs)#, failed_designs

types = ['wing', 'horizontal', 'vertical']
size_str = np.arange(0.03, 0.06, 0.01)
n_str_top = np.arange(2, 8, 1)
n_str_bot = np.arange(1, 6, 1)
t_skin = np.arange(0.002, 0.003, 0.0005)
for type in types:
    if type == 'vertical':
        n = [n_ult_pos, -n_ult_pos, n_ult_pos, n_ult_pos]# for the vertical tail
    else:
        n = [n_ult_pos, n_ult_neg, n_ult_flaps, n_ult_ail] # for the wing
    print(n)
    materials = [AL7055, AL2099]
    start = time()
    print(type, iterate_wing_box(t_skin, n_str_top, n_str_bot, size_str, materials, n, type))
    end = time()
    print(f'Time spend: {end-start}s')
    #print(iterate_wing_box(t_skin, n_str, size_str, AL, n, type))
    #print(iterate_wing_box(t_skin, n_str, size_str, Ti, n, type))
    #print(iterate_wing_box(t_skin, n_str, size_str, Glare, n, type))