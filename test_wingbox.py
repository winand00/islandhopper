from wing_box import wingbox, stringer, skin, crosssection, Material
from design_loads import design_loads, tail_load_elevator
import wingbox_inputs as wb
import unittest

n_max_pos, n_max_neg, _ = design_loads()
n_ult_pos, n_ult_neg = 1.5 * n_max_pos, 1.5 * n_max_neg

t_skin = 0.006
stringers_top = 8
stringers_bot = 2
size_str = 0.03
n_str = 3

# Aluminum 7040
density = 2820  # kg/m3, density of aluminium
E = 69 * 10 ** 9
G = 26.4 * 10 ** 9
sigma_y = 450 * 10 ** 6
poisson = 0.33
AL7040 = Material(density, E, G, sigma_y, poisson)


def make_wingbox(t_skin, n_str, str_size, material, n, type):
    L_D = 10
    w_ac = 84516
    b_wing = 2
    b_vert = 2
    b_hor = 2
    if type == 'wing':
        b = b_wing
        h_box = 1
        w_box = 1
        taper = 1

    if type == 'vertical':
        b = b_vert
        h_box = 1
        w_box = 1
        taper = 1

    if type == 'horizontal':
        b = b_hor
        h_box = wb.h_hor
        w_box = wb.w_hor
        taper = wb.taper_hor

    # wing parameters
    w_L_D = L_D
    w_ly_e = 0.3
    w_lz_e = 0.5
    w_w_engine = 500#(wb.w_engine + wb.w_radiator + wb.w_prop) * 9.81
    powersetting = wb.powersetting
    w_T_engine = 500#w_ac / w_L_D / 2 / powersetting
    w_ly_fc = 0.3
    w_w_fc = 0#wb.w_fc
    w_ly_bat = 0.2
    w_w_bat = 0#wb.w_bat
    w_ly_hld = 0.3
    w_lx_hld = 0.5
    w_F_hld = 200#wb.F_hld
    w_ly_el = 0.8
    w_lx_el = 0.5
    w_F_el =  100#wb.F_a
    w_w_wing = 1000
    w_Mh = 0
    w_lz_h = 0
    w_F_h = 0

    # horizontal tail parameters
    h_L_D = L_D
    h_ly_e = 0
    h_lz_e = 0
    h_w_engine = 0
    h_T_engine = 0
    h_ly_fc = 0
    h_w_fc = 0
    h_ly_bat = 0
    h_w_bat = 0
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

    # vertical tail parameters
    v_L_D = L_D
    v_ly_e = 0
    v_lz_e = 0
    v_w_engine = 0
    v_T_engine = 0
    v_ly_fc = 0
    v_w_fc = 0
    v_ly_bat = 0
    v_w_bat = 0
    v_ly_hld = 0
    v_lx_hld = 0
    v_F_hld = 0  # 15000
    v_ly_el = wb.ly_r
    v_lx_el = wb.lx_r
    v_F_el = 0#wb.F_r
    v_w_wing = 0#wb.L_ver
    v_Mh = 1000#h_F_el * h_lx_el
    v_lz_h = 0.5#wb.lz_hor
    v_F_h = 500# h_w_wing * b_hor

    # skin(height, width, x_coordinate, z_coordinate)
    # coordinates are the bottom left point of the skin
    skin_top = skin(t_skin, w_box, 0, h_box)
    skin_bottom = skin(t_skin, w_box, 0, 0)
    skin_left = skin(h_box, t_skin, 0, 0)
    skin_right = skin(h_box, t_skin, w_box, 0)
    skins = [skin_top, skin_bottom, skin_left, skin_right]

    l_w = b / 2  # length of the wingbox

    # Creating stringerss

    # Y_end position, number of stringers
    y_stringers_stop_top = [(0.2 * l_w, n_str),
                            (0.4 * l_w, n_str),
                            (0.6 * l_w, n_str),
                            (0.8 * l_w, n_str),
                            (1.0 * l_w, n_str)]

    # Y_end position, number of stringers
    y_stringers_stop_bot = [(0.2 * l_w, 0),
                            (0.4 * l_w, 0),
                            (0.6 * l_w, n_str),
                            (0.8 * l_w, 0),
                            (1.0 * l_w, n_str)]
    stringer_width = str_size
    stringer_height = str_size
    stringer_t = 0.005
    stringer_list = {}
    top_stringer_list = []
    bottom_stringer_list = []
    for i in y_stringers_stop_top:
        for j in range(i[1]):
            top_stringer_list.append(stringer(stringer_width, stringer_height, stringer_t, i[0]))
    stringer_list['top'] = top_stringer_list
    for i in y_stringers_stop_bot:
        for j in range(i[1]):
            bottom_stringer_list.append(stringer(stringer_width, stringer_height, stringer_t, i[0]))
    stringer_list['bottom'] = bottom_stringer_list

    root_crosssection = crosssection(stringer_list, skins)

    density_AL = material.density  # kg/m3, density of aluminium
    E = material.E
    G = material.G
    sigma_y = material.sigma_y
    poisson = material.poisson

    if type == 'wing':
        ly_e = w_ly_e
        lz_e = w_lz_e
        w_engine = w_w_engine
        T_engine = w_T_engine
        ly_fc = w_ly_fc
        w_fc = w_w_fc
        ly_bat = w_ly_bat
        w_bat = w_w_bat
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

    return wingbox(stringer_list, root_crosssection, l_w, taper, density_AL, E, G, sigma_y, poisson, ly_e, w_wing,
                   w_engine, L_D,
                   lz_e, ly_hld, lx_hld, ly_el, lx_el, T_engine, F_hld, F_el, Mh, F_h, lz_h, ly_fc, w_fc, ly_bat, w_bat,
                   n, type)

type = 'wing'
wingbox_skins = make_wingbox(0.005, 0, 0.05, AL7040, n_ult_pos, type)
wingbox_stringers = make_wingbox(0.005, 2, 0.1, AL7040, n_ult_pos, type)
wingbox_vertical = make_wingbox(0.005, 0, 0.1, AL7040, n_ult_pos, 'vertical')

class WingboxTest(unittest.TestCase):
    def test_skin_area(self):
        area = wingbox_skins .local_crosssection(0).local_skins[0].area
        self.assertAlmostEqual(area, 0.005, 4)

    def test_stringer_count(self):
        count = len(wingbox_stringers.local_crosssection(0).local_stringers['top'])
        count += len(wingbox_stringers.local_crosssection(0).local_stringers['bottom'])
        self.assertAlmostEqual(count, 14, 4)

    def test_stringer_area(self):
        area = wingbox_stringers.local_crosssection(0).local_stringers['top'][0].area
        self.assertAlmostEqual(area, 9.75*10**-4, 4)

    def test_cross_section_skin_area(self):
        area = wingbox_skins.local_crosssection(0.5).skin_area()
        self.assertAlmostEqual(area, 0.020, 4)

    def test_cross_section_stringer_area(self):
        area = wingbox_stringers.local_crosssection(0).area
        self.assertAlmostEqual(area, 0.03365, 4)

    def test_skin_inertia(self):
        Ixx = wingbox_skins .local_crosssection(0).local_skins[0].I_xx()
        Izz = wingbox_skins .local_crosssection(0).local_skins[0].I_xx()
        self.assertAlmostEqual(Ixx, 1.041667*10 **-8, 3)
        self.assertAlmostEqual(Izz, 4.1667 *10 **-4, 3)

    def test_centroid(self):
        centroid_z = wingbox_skins .local_crosssection(0).centroid_z()
        centroid_x = wingbox_skins .local_crosssection(0).centroid_x()
        self.assertAlmostEqual(centroid_z, 0.5, 2)
        self.assertAlmostEqual(centroid_x, 0.5, 2)

    def test_centroid_stringers(self):
        centroid_z = wingbox_stringers.local_crosssection(0).centroid_z()
        centroid_x = wingbox_stringers.local_crosssection(0).centroid_x()
        self.assertAlmostEqual(centroid_z, 0.5869, 2)
        self.assertAlmostEqual(centroid_x, 0.5, 2)

    def test_crosssection_inertia(self):
        I_zz = wingbox_skins .local_crosssection(0).I_zz
        I_xx = wingbox_skins .local_crosssection(0).I_xx
        self.assertAlmostEqual(I_zz, 3.33*10**-3, 5)
        self.assertAlmostEqual(I_xx, 3.33*10**-3, 5)

    def test_crosssection_inertia_stringers(self):
        I_zz = wingbox_stringers .local_crosssection(0).I_zz
        I_xx = wingbox_stringers .local_crosssection(0).I_xx
        self.assertAlmostEqual(I_zz, 4.154 *10**-3, 4)
        self.assertAlmostEqual(I_xx, 6.491*10**-3, 5)

    def test_weight(self):
        w = wingbox_skins .weight
        self.assertAlmostEqual(w, 56.4, 3)

    def test_moment_x(self):
        Mx = wingbox_skins.momentx(1)
        self.assertAlmostEqual(Mx, 0, 3)
        Mx = wingbox_skins.momentx(0)
        self.assertAlmostEqual(Mx, 73.45, 2)


    def test_moment_z(self):
        Mz = wingbox_skins.momentz(1)
        self.assertAlmostEqual(Mz, 0, 3)
        Mz = wingbox_skins.momentz(0)
        self.assertAlmostEqual(Mz, 100, 3)

    def test_shear_x(self):
        Vx = wingbox_skins.shearx(1)
        self.assertAlmostEqual(Vx, 0, 3)
        Vx = wingbox_skins.shearx(0)
        self.assertAlmostEqual(Vx, -400, 3)

    def test_shear_z(self):
        Vz = wingbox_skins.shearz(1)
        self.assertAlmostEqual(Vz, 0, 3)
        Vz = wingbox_skins.shearz(0)
        self.assertAlmostEqual(Vz, -53.095, 3)

    def test_displacement_x(self):
        vx = wingbox_skins.total_displacement_x(0)
        self.assertAlmostEqual(vx, 0, 3)
        vx = wingbox_skins.total_displacement_x(1)
        #self.assertAlmostEqual(vx, 5.44*10**-7, 10)

    def test_displacement_z(self):
        vz = wingbox_skins.displacementz(0)
        self.assertAlmostEqual(vz, 0, 3)
        vz = wingbox_skins.displacementz(1)
        self.assertAlmostEqual(vz, 1.549*10**-7, 9)

    def test_displacement_x(self):
        vx = wingbox_skins.displacementx(0)
        self.assertAlmostEqual(vx, 0, 3)
        vx = wingbox_skins.displacementx(1)
        self.assertAlmostEqual(vx, -3.37*10**-8, 9)

    def test_torsion(self):
        T = wingbox_skins.Torsiony(1)
        self.assertAlmostEqual(T, 0, 3)
        T = wingbox_skins.Torsiony(0)
        self.assertAlmostEqual(T, -100, 3)

    def test_torsional_constant(self):
        J = wingbox_stringers.local_crosssection(0).J
        self.assertAlmostEqual(J, 5*10**-3, 15)

    def test_twist(self):
        twist = wingbox_skins.twist(0)
        self.assertAlmostEqual(twist, 0, 10)
        twist = wingbox_skins.twist(1)
        self.assertAlmostEqual(twist, -3.7879*10**-8, 10)

    def test_vertical_moment_horizontal(self):
        Mz = wingbox_vertical.momentz(1)
        self.assertAlmostEqual(Mz, 0, 10)
        Mz = wingbox_vertical.momentz(0)
        self.assertAlmostEqual(Mz, -1000, 10)

    def test_sigma_F(self):
        sigma_F = wingbox_vertical.sigma_F(0)
        self.assertAlmostEqual(sigma_F, 25000, 10)
        sigma_F = wingbox_vertical.sigma_F(1)
        self.assertAlmostEqual(sigma_F, 0, 10)

    def test_vertical_force_horizontal(self):
        sigmaz = wingbox_vertical.max_bending_stress(0) /(10**5)
        self.assertAlmostEqual(sigmaz, 1.75, 2)

        











if __name__ == "__main__":
    unittest.main()