import matplotlib.pyplot as plt
import numpy as np


class crosssection:
    def __init__(self, skins):
        self.skins = skins
        self.area = self.area()
        self.x_centroid = self.centroid_x()
        self.z_centroid = self.centroid_z()
        self.I_xx = self.I_xx()
        self.I_zz = self.I_zz()

    def area(self):
        area = 0
        for skin in self.skins:
            area += skin.area
        return area

    def centroid_x(self):
        area_distance = 0
        for skin in self.skins:
            area_distance += (skin.x + skin.b/2) * skin.area
        x_centroid = area_distance/self.area
        return x_centroid

    def centroid_z(self):
        area_distance = 0
        for skin in self.skins:
            area_distance += (skin.z + skin.h/2) * skin.area
        z_centroid = area_distance/self.area
        return z_centroid

    def I_xx(self):
        I_xx = 0
        for skin in self.skins:
            I_xx += skin.I_xx

            print((skin.z + skin.h/2 - self.z_centroid))
            I_xx += skin.area * (skin.z + skin.h/2 - self.z_centroid) ** 2
        return I_xx

    def I_zz(self):
        I_zz = 0
        for skin in self.skins:
            I_zz += skin.I_zz
            I_zz += skin.area * (skin.x + skin.b/2 -self.x_centroid) ** 2
        return I_zz

    def plot(self):
        for skin in self.skins:
            if skin.b > skin.h:
                t = skin.h *1500
                x = [skin.x, skin.x + skin.b]
                y = [skin.z, skin.z]
            else:
                t = skin.b *1500
                x = [skin.x, skin.x]
                y = [skin.z, skin.z + skin.h]
            plt.plot(x, y, 'black', linewidth= t)
        plt.plot(self.x_centroid, self.z_centroid, 'X')
        plt.show()





class skin:
    ""
    def __init__(self, h, b, x, z):
        self.h = h #Height of the skin
        self.b = b #Width of the skin
        self.x = x
        self.z = z
        self.area = h * b
        self.I_xx = self.I_xx()
        self.I_zz = self.I_zz()

    def I_xx(self):
        return self.b*self.h**3/12

    def I_zz(self):
        return self.h*self.b**3/12

class wingbox:
    def __init__(self, cross_section, length, material_density):
        self.cross_section = cross_section
        self.density = material_density
        self.length = length
        self.weight = self.cross_section.area * self.length * self.density
    def moment(self, l1, l2, w_wing, w_engine, y):
        I = self.cross_section.I_xx
        Ra = -w_wing * l2 + w_engine
        M0 = -(w_engine * l1) + (w_wing * (l2 ** 2) / 2)
        M = -(Ra*y+M0+(w_wing/2)*y**2-(w_engine*Macaulay(y,l1,1)))
        return M
    def graphmoment(self, l1, l2, w_wing, w_engine):
        x = np.arange(0, l2, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.moment(l1, l2, w_wing, w_engine, x[i]))
        plt.plot(x, lst)
        return lst
    def bendingstress(self):
    def displacement(self, E, l1, l2, w_wing, w_engine, y):
        I = self.cross_section.I_xx

        Ra = -w_wing * l2 + w_engine
        M0 = -(w_engine * l1) + (w_wing * (l2 ** 2) / 2)
        v = -((-1 / (E * I)) * ((1 / 6 * Ra * (y ** 3)) + (1 / 2 * M0 * (y ** 2)) + (w_wing / 24 * (y ** 4)) - (w_engine / 6 * Macaulay(y, l1, 3))))
        # v = (-1 / (E * I)) * ((1 / 6 * Ra * y ** 3)  - (w_engine / 6 * Macaulay(y, l1, 3)))
        return v

    def graphdisplacement(self, E, l1, l2, w_wing, w_engine):
        x = np.arange(0, l2, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.displacement(E, l1, l2, w_wing, w_engine, x[i]))
        plt.plot(x, lst)
        return lst


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
crosssection = crosssection(skins)
crosssection.plot()
plt.show()


density_AL = 2712 #kg/m3, density of aluminium
wingbox = wingbox(crosssection, 2, density_AL)
print(wingbox.weight)


# Test values for deflection
l2=2
E=70 * 10 **9

wingbox.graphdisplacement(E,1,l2,1000,1000)
#wingbox.graphmoment(1,l2,1000,1000)
#graphdisplacement(E,8.127179999999999e-05,9,l2,1000,1000)


plt.show()
