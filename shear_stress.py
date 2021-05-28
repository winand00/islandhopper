import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
import numpy as np

Vz = 1
Vx = 1
Izz = 1
Ixx = 1
t = 1
h = 1
w =1

class Stress:
    def __init__(self, h, w, Izz, Ixx, t, Vx, Vz, Mx, Mz, T, x_centroid, z_centroid):
        self.h = h
        self.w = w
        self.Izz = Izz
        self.Ixx = Ixx
        self.t = t
        self.Vx = -Vx
        self.Vz = -Vz
        self.Mx = Mx
        self.Mz = Mz
        self.T = T
        self.x_centroid = x_centroid
        self.z_centroid = z_centroid

    def von_mises(self, x ,z):
        return np.sqrt(self.bendingstress_xz(x, z) ** 2 + 3 * self.shear_total(x, z) **2)
    
    def max_von_mises(self):
        stresses = []
        x, z = self.get_xz(10)
        for i in range(len(x)):
            stresses.append(self.von_mises(x[i],z[i]))
        return max(stresses)

    def bendingstress_x(self,x, z):
        return self.Mx*(z-self.z_centroid)/self.Ixx

    def bendingstress_z(self, x, z):
        return -self.Mz * (x - self.x_centroid) / self.Izz

    def bendingstress_xz(self, x, z):
        return self.Mx*(z-self.z_centroid)/self.Ixx - self.Mz * (x - self.x_centroid) / self.Izz


    def shear_torque(self, x, z):
        A = self.h * self.w
        return self.T/(2*self.t*A)

    def shearstress_x(self, x, z):
        if x == 0: #Left skin
            if x < self.h / 2:
                return -self.q1_x(self.h/2 - z)
            else:
                return self.q1_x(z-self.h/2)
        if z == 0: #bottom skin
            if x < self.w/2:
                return -self.q2_x(x)
            else:
                return -self.q2_x(x)
        if x >= self.w: #right skins
            if x < self.h / 2:
                return -self.q1_x(self.h/2-z)
            else:
                return self.q3_x(self.h -z)
        if z >= self.h: #top skin
            if x < self.w / 2:
                return self.q2_x(x)
            else:
                return self.q2_x(x)

    def shearstress_z(self, x, z):
        if x == 0: #Left skin
            if x < self.h / 2:
                return -self.q2_z(z)
            else:
                return -self.q2_z(z)
        if z == 0: #bottom skin
            if x < self.w/2:
                return -self.q1_z(self.w/2-x)
            else:
                return self.q3_z(self.w-x)
        if x >= self.w: #right skins
            if x < self.h / 2:
                return self.q2_z(self.h-z)
            else:
                return self.q2_z(self.h-z)
        if z >= self.h: #top skin
            if x < self.w / 2:
                return -self.q3_z(x)
            else:
                return self.q1_z(x-self.w/2)

    def shearstress_xz(self, x, z):
        return self.shearstress_x(x, z) + self.shearstress_z(x, z)

    def shear_total(self, x, z):
        return self.shearstress_xz(x, z) + self.shear_torque(x, z)

    def q1_z(self, s):
        return -self.Vz/self.Ixx * self.t * s * self.h/2

    def q2_z(self, s):
        return self.q1_z(self.w/2) - self.Vz/self.Ixx * self.t * (self.h/2 * s - 0.5 * s**2)

    def q3_z(self, s):
        return self.q2_z(self.h) + self.Vz / self.Ixx * self.t * s * self.h / 2


    def q1_x(self, s):
        return self.Vx/self.Izz * self.t * s * self.w/2

    def q2_x(self, s):
        return self.q1_x(self.h/2) + self.Vx/self.Izz * self.t * (self.w/2 * s - 0.5 * s**2)

    def q3_x(self, s):
        return self.q2_x(self.w) - self.Vx / self.Izz * self.t * s * self.w / 2

    def plot_von_mises(self, y):
        self.plot(y, self.von_mises)

    def plot_shear_x(self, y):
        self.plot(y, self.shearstress_x)

    def plot_shear_z(self, y):
        self.plot(y, self.shearstress_z)

    def plot_shear_xz(self, y):
        self.plot(y, self.shearstress_xz)

    def plot_shear_T(self, y):
        self.plot(y, self.shear_torque)

    def plot_bending_stress_x(self, y):
        self.plot(y, self.bendingstress_x)

    def plot_bending_stress_z(self, y):
        self.plot(y, self.bendingstress_z)

    def plot_bending_stress_xz(self, y):
        self.plot(y, self.bendingstress_xz)
        
    def get_xz(self, N):
        x1 = np.arange(0, self.w + self.w /N, self.w /N)
        y1 = np.ones(x1.shape) * self.h

        #Bottom skin
        x2 = np.arange(0, self.w + self.w /N, self.w /N)
        y2 = np.zeros(x2.shape)

        #Right skin
        y3 = np.arange(0, self.h + self.h/N, self.h/N)
        x3 = np.ones(y3.shape) * self.w

        #Left skin
        y4 = np.arange(0, self.h + self.h/N, self.h/N)
        x4 = np.zeros(y4.shape)

        x = np.hstack((x1,np.flip(x3,0),np.flip(x2,0),x4))
        ys = np.hstack((y1,np.flip(y3,0),np.flip(y2,0),y4))
        return x, ys
        



    def plot(self, y, function):
        N = 50
        print(self.h)
        print(self.w)
        # Here are many sets of y to plot vs. x

        x, ys = self.get_xz(N)

        
        # We need to set the plot limits, they will not autoscale
        fig, ax = plt.subplots()
        ax.set_xlim(np.min(x) - 0.2*self.w, np.max(x) +0.2*self.w)
        ax.set_ylim(np.min(ys)-0.2*self.h, np.max(ys)+0.2*self.h)

        z = np.array([function(x[i], ys[i]) for i in range(len(ys))])
        points_for_stack = np.array([x, ys]).T.reshape(-1, 1, 2)
        segments_for_coloring = np.concatenate([points_for_stack[:-1], points_for_stack[1:]], axis=1)
        line_segments = LineCollection(segments_for_coloring,  cmap=plt.get_cmap('jet'), norm=plt.Normalize(z.min(), z.max()))
        line_segments.set_array(z)
        ax.add_collection(line_segments)
        axcb = fig.colorbar(line_segments)
        axcb.set_label('Shear stress')
        ax.set_title('Shear stress over the wing box cross section')
        plt.gca().add_collection(line_segments)  # This allows interactive changing of the colormap.
        
        plt.show()

#stress = Stress(h, w, Izz, Ixx, t, Vx, Vz, Y)
#stress.plot_shear_z(1)