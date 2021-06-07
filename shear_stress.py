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
    def __init__(self, h, w, Izz, Ixx, t, Vx, Vz, Mx, Mz, T, x_centroid, z_centroid, sigma_F):
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
        self.sigma_F = sigma_F

    def von_mises(self, x ,z):
        return np.sqrt(self.bendingstress_xz(x, z) ** 2 + 3 * self.shear_total(x, z) **2)
    
    def max_von_mises(self):
        stresses = []
        x, z = self.get_xz(10)
        for i in range(len(x)):
            stresses.append(self.von_mises(x[i],z[i]))
        return max(stresses)

    def bendingstress_x(self,x, z):
        sigma = self.Mx*(z-self.z_centroid)/self.Ixx
        sigma += self.sigma_F
        return sigma

    def bendingstress_z(self, x, z):
        sigma = -self.Mz * (x - self.x_centroid) / self.Izz
        sigma += self.sigma_F
        return sigma

    def bendingstress_xz(self, x, z):
        sigma = self.Mx*(z-self.z_centroid)/self.Ixx - self.Mz * (x - self.x_centroid) / self.Izz
        sigma += self.sigma_F
        return sigma
        
    def max_bending_xz(self):
        stresses = []
        x, z = self.get_xz(10)
        for i in range(len(x)):
            stresses.append(self.bendingstress_xz(x[i],z[i]))
        return max(stresses)
        
    def shear_torque(self, x, z):
        A = self.h * self.w
        return self.T/(2*self.t*A)

    def shearstress_x(self, x, z):
        if x == 0: #Left skin
            if x < self.h / 2:
                return -self.q1_x(self.h/2 - z)/self.t
            else:
                return self.q1_x(z-self.h/2)/self.t
        if z == 0: #bottom skin
            if x < self.w/2:
                return -self.q2_x(x)/self.t
            else:
                return -self.q2_x(x)/self.t
        if x >= self.w: #right skins
            if x < self.h / 2:
                return -self.q1_x(self.h/2-z)/self.t
            else:
                return self.q3_x(self.h -z)/self.t
        if z >= self.h: #top skin
            if x < self.w / 2:
                return self.q2_x(x)/self.t
            else:
                return self.q2_x(x)/self.t

    def shearstress_z(self, x, z):
        if x == 0: #Left skin
            if x < self.h / 2:
                return -self.q2_z(z)/self.t
            else:
                return -self.q2_z(z)/self.t
        if z == 0: #bottom skin
            if x < self.w/2:
                return -self.q1_z(self.w/2-x)/self.t
            else:
                return self.q3_z(self.w-x)/self.t
        if x >= self.w: #right skins
            if x < self.h / 2:
                return self.q2_z(self.h-z)/self.t
            else:
                return self.q2_z(self.h-z)/self.t
        if z >= self.h: #top skin
            if x < self.w / 2:
                return -self.q3_z(x)/self.t
            else:
                return self.q1_z(x-self.w/2)/self.t

    def shearstress_xz(self, x, z):
        return self.shearstress_x(x, z) + self.shearstress_z(x, z)
    

    def shear_total(self, x, z):
        return self.shearstress_xz(x, z) + self.shear_torque(x, z)
       
    def max_shear_total(self):
        stresses = []
        x, z = self.get_xz(10)
        for i in range(len(x)):
            stresses.append(abs(self.shear_total(x[i],z[i])))
        return max(stresses)
        

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

    def plot_von_mises(self, fig, ax, title):
        self.plot(self.von_mises, fig, ax, title)

    def plot_shear_x(self, fig, ax, title):
        self.plot(self.shearstress_x, fig, ax, title)

    def plot_shear_z(self, fig, ax, title):
        self.plot(self.shearstress_z, fig, ax, title)

    def plot_shear_xz(self, fig, ax, title):
        self.plot(self.shearstress_xz, fig, ax, title)

    def plot_shear_T(self, fig, ax, title):
        self.plot(self.shear_torque, fig, ax, title)

    def plot_shear_total(self, fig, ax, title):
        self.plot(self.shear_total, fig, ax, title)

    def plot_bending_stress_x(self, fig, ax, title):
        self.plot(self.bendingstress_x, fig, ax, title)

    def plot_bending_stress_z(self, fig, ax, title):
        self.plot(self.bendingstress_z, fig, ax, title)

    def plot_bending_stress_xz(self, fig, ax, title):
        self.plot(self.bendingstress_xz, fig, ax, title)
        
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
        



    def plot(self, function, fig, ax, title):
        N = 50
        # Here are many sets of y to plot vs. x

        x, ys = self.get_xz(N)

        
        # We need to set the plot limits, they will not autoscale

        ax.set_xlim(np.min(x) - 0.2*self.w, np.max(x) +0.2*self.w)
        ax.set_ylim(np.min(ys)-0.2*self.h, np.max(ys)+0.2*self.h)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("z [m]")
        z = np.array([function(x[i], ys[i]) for i in range(len(ys))])
        points_for_stack = np.array([x, ys]).T.reshape(-1, 1, 2)
        segments_for_coloring = np.concatenate([points_for_stack[:-1], points_for_stack[1:]], axis=1)
        line_segments = LineCollection(segments_for_coloring,  cmap=plt.get_cmap('jet'), norm=plt.Normalize(z.min(), z.max()))
        line_segments.set_array(z)
        ax.set_title(title)
        axcb = fig.colorbar(line_segments, ax = ax)
        axcb.ax.tick_params(labelsize=12)
        #axcb.set_label('Shear stress')
        bp = ax.add_collection(line_segments)


        #bp = plt.gca().add_collection(line_segments)  # This allows interactive changing of the colormap.
        return bp


#stress = Stress(h, w, Izz, Ixx, t, Vx, Vz, Y)
#stress.plot_shear_z(1)