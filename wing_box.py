import matplotlib.pyplot as plt

class wingbox:
    def __init__(self):
        pass



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




skin_top = skin(0.003, 0.3, 0, 0.1)
skin_bottom = skin(0.003, 0.3, 0 ,0)
skin_left = skin(0.1, 0.003, 0, 0)
skin_right = skin(0.1, 0.003, 0.3, 0)
skins = [skin_top, skin_bottom, skin_left, skin_right]
crosssection = crosssection(skins)
crosssection.plot()
print(crosssection.I_zz)

density_AL = 2712 #kg/m3, density of aluminium
wingbox = wingbox(crosssection, 10, density_AL)
print(wingbox.weight)