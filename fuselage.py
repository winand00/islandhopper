from wing_box import Macaulay
import numpy as np
import matplotlib.pyplot as plt

class Stringer:
    def __init__(self, t, w, h, angle=0, D=0):
        self.t = t
        self.w = w
        self.h = h
        self.y = np.sin(angle) * D/2
        self.z = np.cos(angle) * D/2
        self.area = self.area()
        self.angle = angle
        self.Iyy = self.Iyy()

    def area(self):
        return (self.w + self.h - self.t) * self.t

    def Iyy(self):
        return self.area * self.z ** 2


class Fuselage:
    def __init__(self, stringer, n_str, length, t, D, weight_dist, x_w, F_w, x_t, F_t):
        self.weight_dist = weight_dist
        self.x_w = x_w
        self.F_w = F_w
        self.x_t = x_t
        self.F_t = F_t
        self.length = length
        self.t = t
        self.D = D
        self.stringers = self.make_stringers(stringer, int(n_str))
        self.Iyy = self.Iyy()
        self.max_Vz = self.max_Vz() # 2.6 * 10 **5
        self.max_My = self.max_My() # 1.4*10**6
        self.max_T = 1 *10**5
        self.J = np.pi*self.t*self.D**3/4

    def make_stringers(self, stringer, n_str):
        t = stringer.t
        w = stringer.w
        h = stringer.h
        stringers = []
        str_angle = np.pi/2 / (n_str/4 +1)
        angle = 0
        for i in range(4):
            for j in range(n_str//4):
                angle += str_angle
                stringers.append(Stringer(t, w, h, angle, self.D))
            angle += str_angle
        return stringers

    def V_z(self, x):
        V = -self.weight_dist * x  + self.F_w * Macaulay(x, self.x_w, 0) - self.F_t * Macaulay(x, self.x_t, 0)
        return V
    
    def max_Vz(self):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.V_z(x[i]))
        return max(lst)
        
    def My(self, x):
        M0 = self.weight_dist*self.x_w**2/2 -self.F_t * (self.x_t-self.x_w) - self.weight_dist*(self.length-self.x_w)**2/2
        M = M0 * Macaulay(x, self.x_w, 0) -self.weight_dist * x ** 2 /2 + self.F_w * Macaulay(x, self.x_w, 1) - self.F_t * Macaulay(x, self.x_t, 1)
        return M
    
    def max_My(self):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.My(x[i]))
        return max(lst)

    def Iyy(self):
        Iyy = 0
        for stri in self.stringers:
            Iyy += stri.Iyy
        Iyy += np.pi*self.t*self.D**3/8
        return Iyy

    def bending_stress(self, theta):
        z = np.cos(theta)
        return self.max_My*z/self.Iyy

    def q(self, theta):
        q = self.max_Vz/self.Iyy*np.sin(theta)*self.D**2/4*self.t
        return q

    def shear_stress(self, theta):
        return self.q(theta)/self.t + self.max_T*self.D/2/self.J

    def von_mises(self, theta):
        return np.sqrt(self.bending_stress(theta) ** 2 + 3 * self.shear_stress(theta) ** 2)

    def max_von_mises(self):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.von_mises(x[i]))
        return max(lst)

    def graph(self, function, label, ax):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(function(x[i]))
        ax.set_title(label)
        bp = ax.plot(x, lst)
        return bp

    def graph_circ(self, function, label, ax):
        x = np.arange(0, 2*np.pi, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(function(x[i]))
        ax.set_title(label)
        bp = ax.plot(x, lst)
        return bp

    def graphs(self):
        fig, axs = plt.subplots(3, 2, figsize=(20,10))
        fig.suptitle("Forces, moments of the fuselage")
        self.graph(self.V_z, 'Shear force in z-direction', axs[0,0])
        self.graph(self.My, 'Moment around the y-axis', axs[1, 0])
        self.graph_circ(self.shear_stress, 'Shear stress along the skin', axs[0, 1])
        self.graph_circ(self.bending_stress, 'Bending stress along the skin', axs[1,1])
        self.graph_circ(self.von_mises, 'Von mises stress along the skin', axs[2,1])
        #self.graph(self.momentx, 'Moment around the x-axis', axs[1, 1])
        #self.graph_displacmentz(self.graph_displacmentz, 'Displacement in z-direction', axs[2, 1])
        #self.graph(self.total_displacement_z, 'Displacement in z-direction', axs[2, 1])
        #self.graph(self.Torsiony, 'Torsional moment around the y-axis', axs[1, 2])
        #self.graph(self.total_twist, 'Twist around the y-axis', axs[2, 2])
        plt.show()
        
def create_fuselage():
    weight_ac = 84516
    n = 2.93 * 1.5
    length = 15
    weight_wing = 1000
    weight_dist = (weight_ac - weight_wing) /length * n
    t = 0.003
    D = 2
    x_t = length - 0.1
    F_t = weight_ac * n / 5
    x_w = 15/2
    F_w = weight_ac * n - weight_wing + F_t
    
    #Material
    
    
    #Make stringer
    t_str = 0.005
    w_str = 0.1
    h_str = 0.1
    stringer = Stringer(t_str, w_str, h_str)
    
    n_str = 8 #Has to be a multiple of 4
    
    fuselage = Fuselage(stringer, n_str, length, t, D, weight_dist, x_w, F_w, x_t, F_t)
    fuselage.graphs()

if __name__ == "__main__":
    fuselage = create_fuselage()
    print(fuselage.max_von_mises())