import matplotlib.pyplot as plt
import numpy as np
from shear_stress import Stress


class crosssection:
    def __init__(self, stringers, skins, taper = None, y = None, l_w = None):
        self.y = y
        self.taper = taper
        self.l_w = l_w
        self.stringers = stringers
        self.root_skins = skins
        self.root_width = max([x.b for x in self.root_skins])
        self.root_height = max([x.h for x in self.root_skins])
        if y is not None:
            self.local_stringers = self.get_local_stringers(y)
            self.local_skins = self.get_local_skins()
        else:
            self.local_skins = self.root_skins
            self.local_stringers = self.stringers
        self.local_width = max([x.b for x in self.local_skins])
        self.local_height = max([x.h for x in self.local_skins])
        self.local_t = min([x.h for x in self.local_skins])
        self.area = self.area()
        self.x_centroid = self.centroid_x()
        self.z_centroid = self.centroid_z()
        self.I_xx = self.I_xx()
        self.I_zz = self.I_zz()
        self.J = self.get_J()

    def get_stringer_spacing(self):
        width = self.local_width
        top_stringers = len(self.local_stringers['top'])
        top_spacing = width/(top_stringers+1)
        bottom_stringers = len(self.local_stringers['bottom'])
        bottom_spacing = width/(bottom_stringers+1)
        return top_spacing, bottom_spacing

    def get_local_stringers(self, y):
        new_stringers = {}
        top_stringers = []
        for stri in self.stringers['top']:
            if stri.is_present(y):
                top_stringers.append(stri)
        bottom_stringers = []
        for stri in self.stringers['bottom']:
            if stri.is_present(y):
                bottom_stringers.append(stri)
        new_stringers['top'] = top_stringers
        new_stringers['bottom'] = bottom_stringers
        return new_stringers

    def get_local_skins(self):
        skins_ = []
        for skin_ in self.root_skins:
            b = skin_.b
            h = skin_.h
            x = skin_.x
            z = skin_.z
            if b > skin_.h:
                b = self.skin_taper_width()
            else:
                h = self.skin_taper_height()
            if x > 0:
                x = self.skin_taper_width()
            if z > 0:
                z = self.skin_taper_height()

            skins_.append(skin(h, b, x, z))
        return skins_

    def skin_taper_width(self):
        c_r = self.root_width
        c_t = c_r * self.taper
        return (c_r - c_t) * (self.l_w - self.y) / self.l_w + c_t

    def skin_taper_height(self):
        h_r = self.root_height
        h_t = h_r * self.taper
        return (h_r - h_t) * (self.l_w - self.y) / self.l_w + h_t

    def area(self):
        area = 0
        for sk in self.local_skins:
            area += sk.area
        for stri in self.local_stringers['top']:
            area += stri.area
        for stri in self.local_stringers['bottom']:
            area += stri.area
        return area

    def centroid_x(self):
        area_distance = 0
        for sk in self.local_skins:
            area_distance += (sk.x + sk.b/2) * sk.area
        top_spacing, bottom_spacing = self.get_stringer_spacing()
        x = 0
        for stri in self.local_stringers['top']:
            x += top_spacing
            area_distance += x * stri.area
        x = 0
        for stri in self.local_stringers['bottom']:
            x += bottom_spacing
            area_distance += x * stri.area
        x_centroid = area_distance/self.area
        return x_centroid

    def centroid_z(self):
        area_distance = 0
        for sk in self.local_skins:
            area_distance += (sk.z + sk.h/2) * sk.area
        for stri in self.local_stringers['top']:
            area_distance += self.local_height * stri.area
        z_centroid = area_distance/self.area
        return z_centroid

    def I_xx(self):
        I_xx = 0
        for skin in self.local_skins:
            I_xx += skin.I_xx()
            I_xx += skin.area * (skin.z + skin.h/2 - self.z_centroid) ** 2
        for stri in self.local_stringers['top']:
            I_xx += stri.area * (self.local_height - self.z_centroid)**2
        for stri in self.local_stringers['bottom']:
            I_xx += stri.area * self.z_centroid**2
        return I_xx

    def I_zz(self):
        I_zz = 0
        for skin in self.local_skins:
            I_zz += skin.I_zz()
            I_zz += skin.area * (skin.x + skin.b/2 - self.x_centroid) ** 2

        top_spacing, bottom_spacing = self.get_stringer_spacing()
        x = 0
        for stri in self.local_stringers['top']:
            x += top_spacing
            I_zz += stri.area * (x-self.x_centroid)**2
        x = 0
        for stri in self.local_stringers['bottom']:
            x += bottom_spacing

            I_zz += stri.area * (x-self.x_centroid)**2
        return I_zz

    def get_J(self):
        A = self.local_width * self.local_height
        return 1/((self.local_width * 2 + self.local_height * 2)/(4 * A * self.local_t))

    def graph_stress(self, y, Vx, Vz, Mx, Mz, T):
        stress = Stress(self.local_height , self.local_width , self.I_zz, self.I_xx, self.local_t,
                        Vx, Vz, Mx, Mz, T, self.x_centroid, self.z_centroid)
        fig, axs = plt.subplots(2, 3, figsize=(20,10))
        fig.suptitle(f'Different stresses over the cross section, at y = {y}')
        stress.plot_bending_stress_x(fig, axs[0, 0], f'Bending stress due to Mx')
        stress.plot_bending_stress_z(fig, axs[1, 0], f'Bending stress due to Mz')
        stress.plot_shear_x(fig, axs[0, 1], f'Shear stress due to Vx')
        stress.plot_shear_z(fig, axs[1, 1], f'Shear stress due to Vz')
        stress.plot_shear_total(fig, axs[0, 2], f'Total shear stress')
        stress.plot_von_mises(fig, axs[1, 2], f'Von mises stress')
        plt.show()

    def plot(self):
        stringer_spacing_top, stringer_spacing_bottom = self.get_stringer_spacing()

        x_coord = 0
        for skin in self.local_skins:
            if skin.b > skin.h:
                t = min(skin.h *750, 5)
                x = [skin.x, skin.x + skin.b]
                y = [skin.z, skin.z]
            else:
                t = min(skin.b *750, 5)
                x = [skin.x, skin.x]
                y = [skin.z, skin.z + skin.h]
            plt.plot(x, y, 'black', linewidth= t)

        for stri in self.local_stringers['top']:
            x_coord += stringer_spacing_top
            t = min(stri.t * 750, 5)
            x = [x_coord, x_coord + stri.b]
            y = [self.local_height - stri.t, self.local_height- stri.t]
            plt.plot(x, y, 'black', linewidth=t)
            x = [x_coord, x_coord]
            y = [self.local_height- stri.t, self.local_height - stri.h- stri.t]
            plt.plot(x, y, 'black', linewidth=t)

        x_coord = 0
        for stri in self.local_stringers['bottom']:
            x_coord += stringer_spacing_bottom
            t = min(stri.t * 750, 5)
            x = [x_coord, x_coord + stri.b]
            y = [stri.t, stri.t ]
            plt.plot(x, y, 'black', linewidth=t)
            x = [x_coord, x_coord]
            y = [stri.t , stri.h + stri.t ]
            plt.plot(x, y, 'black', linewidth=t)

        plt.plot(self.centroid_x(), self.centroid_z(), 'X')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()


class skin:
    def __init__(self, h, b, x, z):
        self.h = h  # Height of the skin
        self.b = b  # Width of the skin
        self.x = x
        self.z = z
        self.area = h * b

    def I_xx(self):
        return self.b*self.h**3/12

    def I_zz(self):
        return self.h*self.b**3/12


class stringer:
    def __init__(self, b, h, t, y_end):
        self.b = b
        self.h = h
        self.t = t
        self.area = self.area()
        self.y_end = y_end

    def is_present(self, y):
        return y <= self.y_end

    def area(self):
        return (self.b + self.h - self.t)*self.t


class wingbox:
    def __init__(self, stringers, cross_section, length, taper, material_density, E, G, ly_e, w_wing, w_engine, L_D,
                 lz_e, ly_hld, lx_hld, ly_el, lx_el, T_engine, F_hld, F_el):
        self.stringers = stringers
        self.cross_section = cross_section
        self.density = material_density
        self.E = E
        self.G = G
        self.length = length
        self.skins = self.cross_section.root_skins
        self.taper = taper
        self.weight = self.get_weight()
        self.g= 9.80665
        self.ly_e = ly_e
        self.w_wing = w_wing
        self.w_engine = w_engine
        self.L_D = L_D
        self.lz_e = lz_e
        self.ly_hld = ly_hld
        self.lx_hld = lx_hld
        self.ly_el = ly_el
        self.lx_el = lx_el
        self.T_engine = T_engine
        self.F_hld = F_hld
        self.F_el = F_el

    def local_crosssection(self, y):
        return crosssection(self.stringers, self.skins, self.taper, y, self.length)

    def get_weight(self):
        weight = self.local_crosssection(self.length / 2).area * self.length * self.density
        for stri in self.stringers['top']:
            weight -= (self.length - stri.y_end) * stri.area * self.density
        for stri in self.stringers['bottom']:
            weight -= (self.length - stri.y_end) * stri.area * self.density
        return weight


    def local_area(self, y):
        cross_section = self.local_crosssection(y)
        return cross_section.area

    def local_Ixx(self, y):
        cross_section = self.local_crosssection(y)
        return cross_section.I_xx

    def local_Izz(self, y):
        cross_section = self.local_crosssection(y)
        return cross_section.I_zz

    def lift(self, y):
        return 1000 - y * 10

    def w_steps(self, nr_steps):
        stepsize = self.length / nr_steps
        steps_starts = np.arange(0, self.length, self.length / nr_steps)
        l_w = []
        for i in range(0, len(steps_starts)):
            step = steps_starts[i]
            l_w.append(self.lift(step + stepsize / 2))
        steps_starts = np.append(steps_starts, self.length)
        return steps_starts, l_w

    def shearz(self, y):
        Ra = -self.w_wing * self.length + self.w_engine + self.local_crosssection(self.length / 2).area * self.density * self.g * self.length
        Vz = (Ra+(self.w_wing)*y-(self.w_engine*Macaulay(y,self.ly_e,0))-(self.local_crosssection(self.length / 2).area * self.density * self.g * self.length)*y)
        return Vz

    def shearx(self, y):
        Ra = self.w_wing*(1/self.L_D) * self.length - self.T_engine
        Vx = (Ra-(self.w_wing*(1/self.L_D))*y+(self.T_engine*Macaulay(y,self.ly_e,0)))
        return Vx

    def momentx(self, y):
        Ra = -self.w_wing * self.length + self.w_engine + self.local_crosssection(self.length/2).area*self.density*self.g*self.length
        M0 = -(self.w_engine * self.ly_e) + (self.w_wing * (self.length ** 2) / 2)-(self.local_crosssection(self.length/2).area*self.density*self.g*self.length**2/2)
        M = (Ra*y+M0+(self.w_wing/2)*y**2-(self.local_crosssection(self.length/2).area*self.density*self.g/2)*y**2-(self.w_engine*Macaulay(y,self.ly_e,1)))
        return M

    def momentz(self, y):
        Ra = self.w_wing*(1/self.L_D) * self.length - self.T_engine
        M0 = (self.T_engine * self.ly_e) - (self.w_wing*(1/self.L_D) * (self.length ** 2) / 2)
        M = (Ra*y+M0-(self.w_wing*(1/self.L_D)/2)*y**2+(self.T_engine*Macaulay(y,self.ly_e,1)))
        return M

    def moment2(self, y):
        ylst, wlst = self.w_steps(3)
        lift=[]
        liftmoment=[]
        for i in range(len(ylst)-1):
            lift.append(wlst[i]*(ylst[i+1]-ylst[i]))
            liftmoment.append(wlst[i]*((ylst[i+1]-ylst[i])**2)/2)
        Ra = self.w_engine-sum(lift)
        M0 = -(self.w_engine * self.ly_e) +sum(liftmoment)
        liftmacauley=[]
        for i in range(len(wlst)-1):
            liftmacauley.append((wlst[i+1]-wlst[i])/2*(Macaulay(y,ylst[i+1],2)))
        M = (M0+Ra*y+(wlst[0]/2)*y**2-(self.w_engine*Macaulay(y,self.ly_e,1)))-sum(liftmacauley)
        return M

    def displacementx(self, y):
        I = self.local_crosssection(y).I_zz
        Ra = self.w_wing*(1/self.L_D) * self.length - self.T_engine
        M0 = (self.T_engine * self.ly_e) - (self.w_wing*(1/self.L_D) * (self.length ** 2) / 2)
        v = -1 / (self.E * I)*(Ra*y+M0-(self.w_wing*(1/self.L_D)/24)*y**4+(self.T_engine*Macaulay(y, self.ly_e, 3)/6))
        return v

    def new_displacementx(self, y, step):
        I = self.local_crosssection(y - step/2).I_zz
        Ra = self.w_wing*(1/self.L_D) * self.length - self.T_engine
        M0 = (self.T_engine * self.ly_e) - (self.w_wing*(1/self.L_D) * (self.length ** 2) / 2)
        v = -1 / (self.E * I)*(Ra*y+M0-(self.w_wing*(1/self.L_D)/24)*y**4+(self.T_engine*Macaulay(y, self.ly_e, 3)/6))
        return v/self.length*step

    def displacementz(self, y, step):
        I = self.local_crosssection(y/2).I_xx
        Ra = -self.w_wing * self.length + self.w_engine + self.local_crosssection(
            self.length / 2).area * self.density * self.g * self.length
        M0 = -(self.w_engine * self.ly_e) + (self.w_wing * (self.length ** 2) / 2) - (
                    self.local_crosssection(self.length / 2).area * self.density * self.g * self.length ** 2 / 2)
        v = -((-1 / (self.E * I)) * (((1 / 6) * Ra * (y ** 3)) + (1 / 2 * M0 * (y ** 2)) + ((self.w_wing / 24) * (y ** 4))
                                    - (self.w_engine / 6 * Macaulay(y, self.ly_e, 3))-(((self.local_crosssection(y / 2).area * self.density * self.g * y) / 24) * (y ** 4))))
        return v

    def new_displacementz(self, y, step):
        I = self.local_crosssection(y-step/2).I_xx
        Ra = -self.w_wing * self.length + self.w_engine + self.local_crosssection(
            self.length / 2).area * self.density * self.g * self.length
        M0 = -(self.w_engine * self.ly_e) + (self.w_wing * (self.length ** 2) / 2) - (
                    self.local_crosssection(self.length / 2).area * self.density * self.g * self.length ** 2 / 2)
        v = -((-1 / (self.E * I)) * (((1 / 6) * Ra * (y ** 3)) + (1 / 2 * M0 * (y ** 2)) + ((self.w_wing / 24) * (y ** 4))
                                    - (self.w_engine / 6 * Macaulay(y, self.ly_e, 3))-(((self.local_crosssection(y / 2).area * self.density * self.g * y) / 24) * (y ** 4))))
        return v/self.length*step

    def graph_displacmentz(self, function, label, ax):
        displacements = []
        total_displacements = []
        step = 0.01
        y = np.arange(0, self.length, step)
        #displacements.append(0)
        for i in range(len(y)-1):
            start_displacement = self.displacementz(y[i], step)
            end_displacement = self.displacementz(y[i+1], step)
            displacements.append(end_displacement-start_displacement)
            total_displacements.append(sum(displacements))
        ax.set_title(label)
        bp = ax.plot(y[0:-1], total_displacements)
        return bp

    def total_displacement_z(self, y):
        total_v = 0
        step = 0.05
        y = np.arange(0, y, step)
        for i in y:
            total_v += self.new_displacementz(i, step)
        return total_v

    def total_displacement_x(self, y):
        total_v = 0
        step = 0.05
        y = np.arange(0, y, step)
        for i in y:
            total_v += self.new_displacementx(i, step)
        return total_v




    def displacement2(self, y):
        I = self.local_crosssection(y).I_xx
        ylst, wlst = self.w_steps(3)
        lift = []
        liftmoment = []
        for i in range(len(ylst)-1):
            lift.append(wlst[i]*(ylst[i+1]-ylst[i]))
            liftmoment.append(wlst[i]*((ylst[i+1]-ylst[i])**2)/2)
        Ra = self.w_engine-sum(lift)
        M0 = -(self.w_engine * self.ly_e) + sum(liftmoment)
        liftmacauley = []
        for i in range(len(wlst)-1):
            liftmacauley.append((wlst[i+1]-wlst[i])/24*(Macaulay(y, ylst[i+1], 4)))
        v = (-1 / (self.E * I)) * ((1/2*M0*(y**2))+(Ra/6*(y**3))-(self.w_engine/6*Macaulay(y, self.ly_e, 3))+(wlst[0]/24*(y**4))-sum(liftmacauley))
        return v

    def Torsiony(self, y):
        Ta = self.T_engine*self.lz_e-self.F_el*self.lx_el-self.F_hld*self.lx_hld
        T = Ta+self.F_hld*self.lx_hld*Macaulay(y, self.ly_hld, 0)-self.T_engine*self.lz_e*Macaulay(y, self.ly_e, 0) +\
            self.F_el*self.lx_el*Macaulay(y, self.ly_el, 0)
        return T

    def twist(self, y, step):
        J = self.local_crosssection(y - step/2).J
        Ta = self.T_engine*self.lz_e-self.F_el*self.lx_el-self.F_hld*self.lx_hld
        theta = 1/(self.G * J) * (Ta * y +self.F_hld*self.lx_hld*Macaulay(y, self.ly_hld, 1)-self.T_engine*self.lz_e*Macaulay(y, self.ly_e, 1) +\
            self.F_el*self.lx_el*Macaulay(y, self.ly_el, 1))
        return theta/self.length*step

    def total_twist(self, y):
        total_v = 0
        step = 0.05
        y = np.arange(0, y, step)
        for i in y:
            total_v += self.twist(i, step)
        return total_v

    def plot_crosssection(self, y):
        return self.local_crosssection(y).plot()

    def graph(self, function, label, ax):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(function(x[i]))
        ax.set_title(label)
        bp = ax.plot(x, lst)
        return bp

    def graphs(self):
        fig, axs = plt.subplots(3, 3, figsize=(20,10))
        fig.suptitle("Forces, moments and displacements of the wingbox")
        self.graph(self.shearx, 'Shear force in x-direction', axs[0,0])
        self.graph(self.momentz, 'Moment around the z-axis', axs[1, 0])
        self.graph(self.total_displacement_x, 'Displacement in x-direction', axs[2, 0])
        self.graph(self.shearz, 'Shear force in z-direction', axs[0,1])
        self.graph(self.momentx, 'Moment around the x-axis', axs[1, 1])
        #self.graph_displacmentz(self.graph_displacmentz, 'Displacement in z-direction', axs[2, 1])
        self.graph(self.total_displacement_z, 'Displacement in z-direction', axs[2, 1])
        self.graph(self.Torsiony, 'Torsional moment around the y-axis', axs[1, 2])
        self.graph(self.total_twist, 'Twist around the y-axis', axs[2, 2])
        plt.show()

    def graph_properties(self):
        fig, axs = plt.subplots(2, 2, figsize=(20,10))
        fig.suptitle("Forces, moments and displacements of the wingbox")
        self.graph(self.local_area, 'Area cross section over span', axs[0,0])

        self.graph(self.local_Ixx, 'Ixx over the span', axs[0,1])
        self.graph(self.local_Izz, 'Izz over the span', axs[1, 1])
        plt.show()

    def graph_stress(self, y):
        Vx = self.shearx(y)
        Vz = self.shearz(y)
        Mx = self.momentx(y)
        Mz = self.momentz(y)
        T = self.Torsiony(y)
        self.local_crosssection(y).graph_stress(y, Vx, Vz, Mx, Mz, T)

    def get_max_stress(self):
        y_max = None
        max_stress = 0
        y = np.arange(0, self.length, 0.2)
        for i in y:
            cross_section = self.local_crosssection(i)
            Vx = self.shearx(i)
            Vz = self.shearz(i)
            Mx = self.momentx(i)
            Mz = self.momentz(i)
            t = cross_section.local_t
            w = cross_section.local_width
            h = cross_section.local_height
            Ixx = cross_section.I_xx
            Izz = cross_section.I_zz
            T =  self.Torsiony(i)
            x_centroid = cross_section.x_centroid
            z_centroid = cross_section.z_centroid
            stress = Stress(h, w, Izz, Ixx, t, Vx, Vz,Mx, Mz, T, x_centroid, z_centroid)
            if stress.max_von_mises() > max_stress:
                max_stress = stress.max_von_mises()
                y_max = i
        return y_max, max_stress
            

def Macaulay(x, x_point, power):
    if (x-x_point) < 0:
        return 0
    elif power == 0:
        return 1
    else:
        return (x-x_point)**power


b = 20
n = 3.8
w_ac = 84516

t_skin = 0.006
h_box = 0.4
w_box = 1.5

# skin(height, width, x_coordinate, z_coordinate)
# coordinates are the bottom left point of the skin
skin_top = skin(t_skin, w_box, 0, h_box)
skin_bottom = skin(t_skin, w_box, 0, 0)
skin_left = skin(h_box, t_skin, 0, 0)
skin_right = skin(h_box, t_skin, w_box, 0)
skins = [skin_top, skin_bottom, skin_left, skin_right]
taper = 0.5  # taper ratio of the wingbox
l_w = b/2  # length of the wingbox
y_e = 1  # location of the engine

# crosssection.plot()
# plt.show()

# Creating stringerss

# Y_end position, number of stringers
y_stringers_stop_top = [(0.2 * l_w, 4),
                        (0.4 * l_w, 4),
                        (0.6 * l_w, 4),
                        (0.8 * l_w, 2),
                        (1.0 * l_w, 2)]

# Y_end position, number of stringers
y_stringers_stop_bot = [(0.2 * l_w, 0),
                        (0.4 * l_w, 0),
                        (0.6 * l_w, 2),
                        (0.8 * l_w, 0),
                        (1.0 * l_w, 2)]

stringer_width = 0.05
stringer_height = 0.05
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
print(bottom_stringer_list)
root_crosssection = crosssection(stringer_list, skins)
density_AL = 2712  # kg/m3, density of aluminium

# Test values for deflection

E = 70 * 10 ** 9
G = 26 * 10 ** 9

# Variables from other departments
ly_e = 0.2 * b/2
lz_e = h_box/2
ly_hld = 0.3 * b/2
lx_hld = w_box/2
ly_el = 0.8 * b/2
lx_el = w_box/2
w_engine = 200 * 9.81
w_wing = w_ac / b * n
T_engine = 1300*1000 / 90
F_el = 0
F_hld = 15000
L_D = 12


wingbox = wingbox(stringer_list, root_crosssection, l_w, taper, density_AL, E, G, ly_e, w_wing, w_engine, L_D,
                  lz_e, ly_hld, lx_hld, ly_el, lx_el, T_engine, F_hld, F_el)

wingbox.plot_crosssection(5)
plt.show()

y_max, max_stress = wingbox.get_max_stress()
print(f'{max_stress/(10**6)} MPa, at y = {y_max} m')
print(wingbox.weight, 'kg')
wingbox.graph_stress(y_max)
wingbox.graphs()
wingbox.graph_properties()
