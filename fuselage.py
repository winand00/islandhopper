from wing_box import Macaulay
import numpy as np


class Fuselage:
    def __init__(self, weight_dist, x_wing):
        self.weight_dist
        self.x_wing

    def shear_force_z(self, x):
        V = -self.weight_dist * x**2 / 2 + F_w * Macaulay(x, self.x_wing, 0) - F_t * Macaulay(x, x_t, 0)
        return V

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
        self.graph(self.shearforce_z, 'Shear force in x-direction', axs[0,0])
        #self.graph(self.momentz, 'Moment around the z-axis', axs[1, 0])
        #self.graph(self.total_displacement_x, 'Displacement in x-direction', axs[2, 0])
        #self.graph(self.shearz, 'Shear force in z-direction', axs[0,1])
        #self.graph(self.momentx, 'Moment around the x-axis', axs[1, 1])
        #self.graph_displacmentz(self.graph_displacmentz, 'Displacement in z-direction', axs[2, 1])
        #self.graph(self.total_displacement_z, 'Displacement in z-direction', axs[2, 1])
        #self.graph(self.Torsiony, 'Torsional moment around the y-axis', axs[1, 2])
        #self.graph(self.total_twist, 'Twist around the y-axis', axs[2, 2])
        plt.show()