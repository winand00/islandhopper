from wing_box import wingbox, stringer, skin, crosssection, Material
from design_loads import design_loads, tail_load_elevator
import wingbox_inputs as wb

n_max_pos, n_max_neg, n_max_flaps = design_loads()
n_ult_pos, n_ult_neg, n_ult_flaps = 1.5*n_max_pos, 1.5*n_max_neg, 1.5 * n_max_flaps
n_ult_ail = n_ult_pos * 2/3
t_skin = 0.006
stringers_top = 8
stringers_bot = 2
size_str = 0.03
n_str = 3


#Aluminum 7040
density = 2820  # kg/m3, density of aluminium
E = 69 * 10 ** 9
G = 26.4 * 10 ** 9
sigma_y = 450 * 10 ** 6
poisson = 0.33
AL7040 = Material(density, E, G, sigma_y, poisson, 'AL7040')

#Aluminum 
density = 2810  # kg/m3, density of aluminium
E = 71.1 * 10 ** 9
G = 26.9 * 10 ** 9
sigma_y = 550 * 10 ** 6
poisson = 0.33
AL = Material(density, E, G, sigma_y, poisson, 'AL')

#Aluminium 7055
density = 2860  # kg/m3, density of aluminium
E = 75 * 10 ** 9
G = 26.9 * 10 ** 9
sigma_y = 614 * 10 ** 6
poisson = 0.33
AL7055 = Material(density, E, G, sigma_y, poisson, 'AL7055')

#Aluminium 2099
density = 2630  # kg/m3, density of aluminium
E = 78 * 10 ** 9
G = 27.5 * 10 ** 9
sigma_y = 455 * 10 ** 6
poisson = 0.33
AL2099 = Material(density, E, G, sigma_y, poisson, 'AL2099')

#Aluminium 6061
density = 2700  # kg/m3, density of aluminium
E = 69 * 10 ** 9
G = 26.2 * 10 ** 9
sigma_y = 455 * 10 ** 6
poisson = 0.33
AL6061 = Material(density, E, G, sigma_y, poisson, 'AL6061')

#Titanium
density = 4800  # kg/m3, density of aluminium
E = 110 * 10 ** 9
G = 90 * 10 ** 9
sigma_y = 1200 * 10 ** 6
poisson = 0.31
Ti = Material(density, E, G, sigma_y, poisson, 'Titanium')

#Glare
density = 2520  # kg/m3, density of aluminium
E = 67.9 * 10 ** 9
G = 26.7 * 10 ** 9
sigma_y = 280 * 10 ** 6
poisson = 0.33
Glare = Material(density, E, G, sigma_y, poisson, 'Glare')




def make_wingbox(t_skin, n_str_top, n_str_bot, str_size, material_skin, material_stringers_top, material_stringers_bottom, n, type):
    if n != n_ult_flaps:
        flaps_on = 0
    else:
        flaps_on = 1
    if n != n_ult_ail:
        ail_on = 0
    else:
        ail_on = 1


    L_D = wb.L_D
    w_ac = 84516
    b_wing = wb.b_wing
    b_vert = wb.b_vert
    b_hor = wb.b_hor
    if type == 'wing':
        b = b_wing
        h_box = wb.h_wing
        w_box = wb.w_wing
        taper = 1#wb.taper_wing

    if type == 'vertical':
        b = b_vert
        h_box = wb.h_vert
        w_box = wb.w_vert
        taper = wb.taper_vert

    if type == 'horizontal':
        b = b_hor
        h_box = wb.h_hor
        w_box = wb.w_hor
        taper = wb.taper_hor

    #wing parameters
    w_L_D = L_D
    w_ly_e = wb.ly_e
    w_lz_e = wb.lz_e
    w_w_engine = (wb.w_engine + wb.w_prop)
    powersetting = wb.powersetting
    w_T_engine = w_ac / w_L_D / 2 / powersetting
    w_ly_fc = wb.ly_fc
    w_w_fc = wb.w_fc
    w_ly_bat = wb.ly_bat
    w_w_bat = wb.w_bat
    w_ly_sys = wb.ly_sys
    w_w_sys = wb.w_sys
    w_ly_hld = wb.ly_hld
    w_lx_hld = wb.lx_hld
    w_F_hld = flaps_on * wb.F_hld
    w_ly_el = wb.ly_a
    w_lx_el = wb.lx_a
    w_F_el = ail_on * wb.F_a
    w_w_wing = w_ac / b * n
    w_Mh = 0
    w_lz_h = 0
    w_F_h = 0

    #horizontal tail parameters
    h_L_D = L_D
    h_ly_e = 0
    h_lz_e = 0
    h_w_engine = 0
    h_T_engine = 0
    h_ly_fc = 0
    h_w_fc = 0
    h_ly_bat = 0
    h_w_bat = 0
    h_ly_sys = 0
    h_w_sys = 0
    h_ly_hld = 0
    h_lx_hld = 0
    h_F_hld = 0
    h_ly_el = wb.ly_el  # y position of the elevator
    h_lx_el = wb.lx_el  # x position of the elevator
    h_F_el = wb.F_el  # elevator force
    h_w_wing = wb.L_hor * n  # Is scaled with the loadfactor
    h_Mh = 0
    h_lz_h = 0
    h_F_h = 0

    #vertical tail parameters
    v_L_D = L_D
    v_ly_e = 0
    v_lz_e = 0
    v_w_engine = 0
    v_T_engine = 0
    v_ly_fc = 0
    v_w_fc = 0
    v_ly_bat = 0
    v_w_bat = 0
    v_ly_sys = 0
    v_w_sys = 0
    v_ly_hld = 0
    v_lx_hld = 0
    v_F_hld = 0  # 15000
    v_ly_el = wb.ly_r
    v_lx_el = wb.lx_r
    v_F_el = wb.F_r
    v_w_wing = wb.L_ver
    v_Mh = h_F_el * h_lx_el
    v_lz_h = wb.lz_hor
    v_F_h = -h_w_wing * b_hor
     

    # skin(height, width, x_coordinate, z_coordinate)
    # coordinates are the bottom left point of the skin
    skin_top = skin(t_skin, w_box, 0, h_box)
    skin_bottom = skin(t_skin, w_box, 0, 0)
    skin_left = skin(h_box, t_skin, 0, 0)
    skin_right = skin(h_box, t_skin, w_box, 0)
    skins = [skin_top, skin_bottom, skin_left, skin_right]

    l_w = b/2  # length of the wingbox

    # Creating stringerss
    #material_stringers_top = AL
    #material_stringers_bottom = AL7040
    # Y_end position, number of stringers
    y_stringers_stop_top = [(0.2 * l_w, n_str_top),
                            (0.4 * l_w, n_str_top),
                            (0.6 * l_w, n_str_top),
                            (0.8 * l_w, n_str_top),
                            (1.0 * l_w, n_str_top)]

    # Y_end position, number of stringers
    y_stringers_stop_bot = [(0.2 * l_w, 0),
                            (0.4 * l_w, 0),
                            (0.6 * l_w, n_str_bot),
                            (0.8 * l_w, 0),
                            (1.0 * l_w, n_str_bot)]
    stringer_width = str_size
    stringer_height = str_size
    stringer_t = 0.005
    stringer_list = {}
    top_stringer_list = []
    bottom_stringer_list = []
    for i in y_stringers_stop_top:
        for j in range(i[1]):
            top_stringer_list.append(stringer(stringer_width, stringer_height, stringer_t, i[0], material_stringers_top))
    stringer_list['top'] = top_stringer_list
    for i in y_stringers_stop_bot:
        for j in range(i[1]):
            bottom_stringer_list.append(stringer(stringer_width, stringer_height, stringer_t, i[0], material_stringers_bottom))
    stringer_list['bottom'] = bottom_stringer_list

    root_crosssection = crosssection(stringer_list, skins)

    density_AL = material_skin.density # kg/m3, density of aluminium
    E = material_skin.E
    G = material_skin.G
    sigma_y = material_skin.sigma_y
    poisson = material_skin.poisson

    if type == 'wing':
        ly_e = w_ly_e
        lz_e = w_lz_e
        w_engine = w_w_engine
        T_engine = w_T_engine
        ly_fc = w_ly_fc
        w_fc = w_w_fc
        ly_bat = w_ly_bat
        w_bat = w_w_bat
        ly_sys = w_ly_sys
        w_sys = w_w_sys
        ly_hld = w_ly_hld
        lx_hld = w_lx_hld
        F_hld = w_F_hld
        ly_el = w_ly_el
        lx_el = w_lx_el
        F_el = w_F_el
        w_wing = w_w_wing
        L_D = w_L_D
        Mh = w_Mh
        lz_h = w_lz_h
        F_h = w_F_h

    if type == 'horizontal':
        ly_e = h_ly_e
        lz_e = h_lz_e
        w_engine = h_w_engine
        T_engine = h_T_engine
        ly_fc = h_ly_fc
        w_fc = h_w_fc
        ly_bat = h_ly_bat
        w_bat = h_w_bat
        ly_sys = h_ly_sys
        w_sys = h_w_sys
        ly_hld = h_ly_hld
        lx_hld = h_lx_hld
        F_hld = h_F_hld
        ly_el = h_ly_el
        lx_el = h_lx_el
        F_el = h_F_el
        w_wing = h_w_wing
        L_D = h_L_D
        Mh = h_Mh
        lz_h = h_lz_h
        F_h = h_F_h

    if type == 'vertical':
        ly_e = v_ly_e
        lz_e = v_lz_e
        w_engine = v_w_engine
        T_engine = v_T_engine
        ly_fc = v_ly_fc
        w_fc = v_w_fc
        ly_bat = v_ly_bat
        w_bat = v_w_bat
        ly_sys = v_ly_sys
        w_sys = v_w_sys
        ly_hld = v_ly_hld
        lx_hld = v_lx_hld
        F_hld = v_F_hld
        ly_el = v_ly_el
        lx_el = v_lx_el
        F_el = v_F_el
        w_wing = v_w_wing
        L_D = v_L_D
        Mh = v_Mh
        lz_h = v_lz_h
        F_h = v_F_h



    return wingbox(stringer_list, root_crosssection, l_w, taper, density_AL, E, G, sigma_y, poisson, ly_e, w_wing, w_engine, L_D,
                      lz_e, ly_hld, lx_hld, ly_el, lx_el, T_engine, F_hld, F_el, Mh, F_h, lz_h, ly_fc, w_fc, ly_bat, w_bat, ly_sys, w_sys, n, type)


if __name__ == "__main__":
    type = 'wing'
    wingbox = make_wingbox(0.004, 1, 1, 0.03, AL7040, AL7040, AL7040, n_ult_pos, type)
    print(wingbox.front_skin_weight())
    wingbox.plot_crosssection(0)
    # plt.show()
    print(wingbox.get_max_stress())
    # wingbox.graph_properties()
    wingbox.graphs()
    print(wingbox.top_stringer_weight(), wingbox.bottom_stringer_weight(), wingbox.skin_weight())
    y_max, max_stress = wingbox.get_max_stress()
    print(f'{max_stress/(10**6)} MPa, at y = {y_max} m')
    print(wingbox.weight, 'kg')
    # wingbox.graph_stress(5)

    #print(wingbox.max_bending_stress(0))
    #print(wingbox.is_buckling())
    #print(wingbox.is_crippling())