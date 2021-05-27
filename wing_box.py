import matplotlib.pyplot as plt
import numpy as np


class crosssection:
    def __init__(self, stringers, skins, taper = None, y = None, l_w = None):
        self.y = y
        self.taper = taper
        self.l_w = l_w
        self.stringers = stringers
        self.root_skins = skins
        self.root_width = max([x.b for x in self.root_skins])
        if y is not None:
            self.local_skins = self.get_local_skins()
        else:
            self.local_skins = self.root_skins

        self.local_width = max([x.b for x in self.local_skins])
        self.local_height = max([x.h for x in self.local_skins])
        self.area = self.area()
        self.x_centroid = self.centroid_x()
        self.z_centroid = self.centroid_z()
        self.I_xx = self.I_xx()
        self.I_zz = self.I_zz()

    def get_stringer_spacing(self):
        width = self.local_width
        top_stringers = len(self.stringers['top'])
        top_spacing = width/(top_stringers+1)
        bottom_stringers = len(self.stringers['bottom'])
        bottom_spacing = width/(bottom_stringers+1)
        return top_spacing, bottom_spacing



    def get_local_skins(self):
        skins_ = []
        for skin_ in self.root_skins:
            b = skin_.b
            x = skin_.x
            if b > skin_.h:
                b = self.skin_taper()
            if x > 0:
                x = self.skin_taper()
            skins_.append(skin(skin_.h, b, x, skin_.z))
        return skins_

    def skin_taper(self):
        c_r = self.root_width
        c_t = c_r * self.taper
        return (c_r - c_t) * (self.l_w - self.y) / self.l_w + c_t


    def area(self):
        area = 0
        for sk in self.local_skins:
            area += sk.area
        for stri in self.stringers['top']:
            area += stri.area
        for stri in self.stringers['bottom']:
            area += stri.area
        return area

    def centroid_x(self):
        area_distance = 0
        for sk in self.local_skins:
            area_distance += (sk.x + sk.b/2) * sk.area
        top_spacing, bottom_spacing = self.get_stringer_spacing()
        x = 0
        for stri in self.stringers['top']:
            x += top_spacing
            area_distance += x * stri.area
        x = 0
        for stri in self.stringers['bottom']:
            x += bottom_spacing
            area_distance += x * stri.area
        x_centroid = area_distance/self.area
        return x_centroid

    def centroid_z(self):
        area_distance = 0
        for sk in self.local_skins:
            area_distance += (sk.z + sk.h/2) * sk.area
        for stri in self.stringers['top']:
            area_distance += self.local_height * stri.area
        z_centroid = area_distance/self.area
        return z_centroid

    def I_xx(self):
        I_xx = 0
        for skin in self.local_skins:
            I_xx += skin.I_xx()
            I_xx += skin.area * (skin.z + skin.h/2 - self.z_centroid) ** 2
        for stri in self.stringers['top']:
            I_xx += stri.area * (self.local_height - self.z_centroid)**2
        for stri in self.stringers['bottom']:
            I_xx += stri.area * self.z_centroid**2
        return I_xx

    def I_zz(self):
        I_zz = 0
        for skin in self.local_skins:
            I_zz += skin.I_zz()
            I_zz += skin.area * (skin.x + skin.b/2 - self.x_centroid) ** 2

        top_spacing, bottom_spacing = self.get_stringer_spacing()
        x = 0
        for stri in self.stringers['top']:
            x += top_spacing
            I_zz += stri.area * (x-self.x_centroid)**2
        x = 0
        for stri in self.stringers['bottom']:
            x += top_spacing
            I_zz += stri.area * (x-self.x_centroid)**2

        return I_zz

    def plot(self):
        stringer_spacing_top, stringer_spacing_bottom = self.get_stringer_spacing()
        x_coord = 0
        for skin in self.local_skins:
            if skin.b > skin.h:
                t = skin.h *1500
                x = [skin.x, skin.x + skin.b]
                y = [skin.z, skin.z]
            else:
                t = skin.b *1500
                x = [skin.x, skin.x]
                y = [skin.z, skin.z + skin.h]
            plt.plot(x, y, 'black', linewidth= t)

        for stri in self.stringers['top']:
            x_coord += stringer_spacing_top
            t = stri.t * 1500
            x = [x_coord, x_coord + stri.b]
            y = [self.local_height - stri.t/2, self.local_height- stri.t/2]
            plt.plot(x, y, 'black', linewidth=t)
            x = [x_coord, x_coord]
            y = [self.local_height- stri.t/2, self.local_height - stri.h- stri.t/2]
            plt.plot(x, y, 'black', linewidth=t)
        x_coord = 0
        for stri in self.stringers['bottom']:
            x_coord += stringer_spacing_bottom
            t = stri.t * 1500
            x = [x_coord, x_coord + stri.b]
            y = [stri.t/2, stri.t/2 ]
            plt.plot(x, y, 'black', linewidth=t)
            x = [x_coord, x_coord]
            y = [stri.t/2 , stri.h + stri.t/2 ]
            plt.plot(x, y, 'black', linewidth=t)

        plt.plot(self.centroid_x(), self.centroid_z(), 'X')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()





class skin:
    ""
    def __init__(self, h, b, x, z):
        self.h = h #Height of the skin
        self.b = b #Width of the skin
        self.x = x
        self.z = z
        self.area = h * b


    def I_xx(self):
        b = self.b
        return self.b*self.h**3/12

    def I_zz(self):
        b = self.b
        return self.h*self.b**3/12

class stringer:
    def __init__(self, b, h, t):
        self.b = b
        self.h = h
        self.t = t
        self.area = self.area()

    def area(self):
        return (self.b + self.h - self.t)*self.t


class wingbox:
    def __init__(self, stringers, cross_section, length, taper, material_density):
        self.stringers = stringers
        self.cross_section = cross_section
        self.density = material_density
        self.length = length
        #self.weight = self.cross_section.area * self.length * self.density
        self.skins = self.cross_section.root_skins
        self.taper = taper
        self.g=9.80665
    def local_crosssection(self, y):
        return crosssection(self.stringers, self.skins, self.taper, y, self.length)

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

    def shearz(self, l1, w_wing, w_engine, y):
        Ra = -w_wing * self.length + w_engine + self.local_crosssection(self.length / 2).area * self.density * self.g * self.length
        P = (Ra+(w_wing)*y-(w_engine*Macaulay(y,l1,0))-(self.local_crosssection(self.length / 2).area * self.density * self.g * self.length)*y)
        return P

    def shearx(self, l1, w_wing, T_engine,L_D, y):
        Ra = w_wing*(1/L_D) * self.length - T_engine
        M = (Ra-(w_wing*(1/L_D))*y+(T_engine*Macaulay(y,l1,0)))
        return M

    def graphshearz(self, l1, w_wing, w_engine):
        x = np.arange(0, self.length, 0.001)
        lst = []
        for i in range(len(x)):
            lst.append(self.shearz(l1, w_wing, w_engine, x[i]))
        plt.plot(x, lst,label='Shear in z direction')
        return lst
    def graphshearx(self, l1, w_wing, w_engine,L_D):
        x = np.arange(0, self.length, 0.001)
        lst = []
        for i in range(len(x)):
            lst.append(self.shearx(l1, w_wing, w_engine, L_D, x[i]))
        plt.plot(x, lst,label='Shear in x direction')
        return lst

    def momentx(self, l1, w_wing, w_engine, y):
        Ra = -w_wing * self.length + w_engine + self.local_crosssection(self.length/2).area*self.density*self.g*self.length
        M0 = -(w_engine * l1) + (w_wing * (self.length ** 2) / 2)-(self.local_crosssection(self.length/2).area*self.density*self.g*self.length**2/2)
        M = (Ra*y+M0+(w_wing/2)*y**2-(self.local_crosssection(self.length/2).area*self.density*self.g/2)*y**2-(w_engine*Macaulay(y,l1,1)))
        return M


    def momentz(self, l1, w_wing, T_engine,L_D, y):
        Ra = w_wing*(1/L_D) * self.length - T_engine
        M0 = (T_engine * l1) - (w_wing*(1/L_D) * (self.length ** 2) / 2)
        M = (Ra*y+M0-(w_wing*(1/L_D)/2)*y**2+(T_engine*Macaulay(y,l1,1)))
        return M

    def moment2(self, l1, w_engine, y):
        ylst, wlst = self.w_steps(3)
        print(ylst)
        print(wlst)
        lift=[]
        liftmoment=[]
        for i in range(len(ylst)-1):
            lift.append(wlst[i]*(ylst[i+1]-ylst[i]))
            liftmoment.append(wlst[i]*((ylst[i+1]-ylst[i])**2)/2)
        Ra = w_engine-sum(lift)
        M0 = -(w_engine * l1) +sum(liftmoment)
        liftmacauley=[]
        for i in range(len(wlst)-1):
            liftmacauley.append((wlst[i+1]-wlst[i])/2*(Macaulay(y,ylst[i+1],2)))
        M = (M0+Ra*y+(wlst[0]/2)*y**2-(w_engine*Macaulay(y,l1,1)))-sum(liftmacauley)
        return M

    def graphmomentx(self, l1, w_wing, w_engine):
        x = np.arange(0, self.length, 0.001)
        lst = []
        for i in range(len(x)):
            lst.append(self.momentx(l1, w_wing, w_engine, x[i]))
        plt.plot(x, lst,label='Moment in x direction')
        return lst
    def graphmomentz(self, l1, w_wing, w_engine,L_D):
        x = np.arange(0, self.length, 0.001)
        lst = []
        for i in range(len(x)):
            lst.append(self.momentz(l1, w_wing, w_engine, L_D, x[i]))
        plt.plot(x, lst,label='Moment in z direction')
        return lst

    def graphmoment2(self, l1, w_engine):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.moment2(l1, w_engine, x[i]))
        plt.plot(x, lst)
        return lst
    def bendingstress(self):
        pass

    def displacementz(self, E, l1, w_wing, w_engine, y):
        I = self.local_crosssection(y).I_xx
        Ra = -w_wing * self.length + w_engine + self.local_crosssection(
            self.length / 2).area * self.density * self.g * self.length
        M0 = -(w_engine * l1) + (w_wing * (self.length ** 2) / 2) - (
                    (self.local_crosssection(self.length / 2).area * self.density * self.g * self.length) ** 2 / 2)
        v = ((-1 / (E * I)) * (((1 / 6) * Ra * (y ** 3)) + (1 / 2 * M0 * (y ** 2)) + ((w_wing / 24) * (y ** 4)) - (w_engine / 6 * Macaulay(y, l1, 3))-(((self.local_crosssection(self.length / 2).area * self.density * self.g * self.length) / 24) * (y ** 4))))
        return v
    def displacement2(self, E, l1, w_engine, y):
        I = self.local_crosssection(y).I_xx
        ylst, wlst = self.w_steps(3)
        lift=[]
        liftmoment=[]
        for i in range(len(ylst)-1):
            lift.append(wlst[i]*(ylst[i+1]-ylst[i]))
            liftmoment.append(wlst[i]*((ylst[i+1]-ylst[i])**2)/2)
        Ra = w_engine-sum(lift)
        M0 = -(w_engine * l1) +sum(liftmoment)
        liftmacauley=[]
        for i in range(len(wlst)-1):
            liftmacauley.append((wlst[i+1]-wlst[i])/24*(Macaulay(y,ylst[i+1],4)))
        v = (-1 / (E * I)) * ((1/2*M0*(y**2))+(Ra/6*(y**3))-(w_engine/6*Macaulay(y,l1,3))+(wlst[0]/24*(y**4))-sum(liftmacauley))
        return v

    def graphdisplacementz(self, E, l1, w_wing, w_engine):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.displacementz(E, l1, w_wing, w_engine, x[i]))
        plt.plot(x, lst,label='Displacement in z direction')
        return lst
    def graphdisplacement2(self, E, l1, w_engine):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.displacement2(E, l1, w_engine, x[i]))
        plt.plot(x, lst)
        return lst
    def plot_crosssection(self, y):
        return self.local_crosssection(y).plot()




def Macaulay(x, x_point, power):
    if (x-x_point) < 0:
        return 0
    elif power == 0:
        return 1
    else:
        return ((x-x_point)**power)



#skin(height, width, x_coordinate, z_coordinate)
#coordinates are the bottom left point of the skin
skin_top = skin(0.003, 0.3, 0, 0.1)
skin_bottom = skin(0.003, 0.3, 0 ,0)
skin_left = skin(0.1, 0.003, 0, 0)
skin_right = skin(0.1, 0.003, 0.3, 0)
skins = [skin_top, skin_bottom, skin_left, skin_right]
taper = 0.5 #taper ratio of the wingbox
l_w = 3 #length of the wingbox
y_e = 1 #location of the engine

#crosssection.plot()
#plt.show()

#Creating stringers
number_top_stringers = 2
number_bottom_stringers = 2
stringer_width = 0.01
stringer_height = 0.01
stringer_t = 0.005
stringer_list = {}
top_stringer_list = []
bottom_stringer_list = []
for i in range(number_top_stringers):
    top_stringer_list.append(stringer(stringer_width, stringer_height, stringer_t))
stringer_list['top'] = top_stringer_list
for i in range(number_bottom_stringers):
    bottom_stringer_list.append(stringer(stringer_width, stringer_height, stringer_t))
stringer_list['bottom'] = bottom_stringer_list

root_crosssection = crosssection(stringer_list, skins)
density_AL = 2712 #kg/m3, density of aluminium
wingbox = wingbox(stringer_list, root_crosssection, l_w, taper, density_AL)
#wingbox.plot_crosssection(0)

#print(wingbox.local_crosssection(0).I_zz)
# Test values for deflection

E=70 * 10 **9
l1=1.5

w_engine=2500
w_wing=300
T_engine=2500
L_D=1

wingbox.graphmomentx(l1,w_wing,w_engine)
wingbox.graphmomentz(l1,w_wing,T_engine,L_D)
wingbox.graphshearz(l1,w_wing,w_engine)
wingbox.graphshearx(l1,w_wing,T_engine,L_D)
wingbox.graphdisplacementz(E,l1,w_wing,w_engine)
plt.legend()
plt.show()

