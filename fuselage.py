from wing_box import Macaulay
import wingbox_inputs as wb
from create_wingbox import AL, AL7040, Glare
import numpy as np
import matplotlib.pyplot as plt

class Stringer:
    def __init__(self, t, w, h, rho, sigma_y, E, poisson, angle=0, D=0):
        self.t = t
        self.w = w
        self.h = h
        self.rho = rho
        self.y = np.sin(angle) * D/2
        self.z = np.cos(angle) * D/2
        self.area = self.area()
        self.angle = angle
        self.Iyy = self.Iyy()
        self.sigma_y = sigma_y
        self.E = E
        self.poisson = poisson

    def area(self):
        return (self.w + self.h - self.t) * self.t

    def Iyy(self):
        return self.area * self.z ** 2
    
    def crippling(self):
        C = 0.425
        n = 0.6
        alpha = 0.8
        h = self.h-self.t
        w = self.w - self.t
        A1 = h * self.t
        A2 = w * self.t
        sigma_cc_1 = self.sigma_y * alpha*(C/self.sigma_y*np.pi**2*self.E/(12*(1-self.poisson**2))*(self.t/h)**2)**(1-n)
        sigma_cc_2 = self.sigma_y * alpha*(C/self.sigma_y*np.pi**2*self.E/(12*(1-self.poisson**2))*(self.t/w)**2)**(1-n)
        sigma_total = (sigma_cc_1*A1+sigma_cc_2*A2)/(A1+A2)
        return sigma_total


class Fuselage:
    def __init__(self, stringer, n_str, length, t, D, weight_dist, x_w, T_t, x_t, F_t, rho, E, sigma_y, poisson, buckling_factor):
        self.weight_dist = weight_dist
        self.rib_spacing = 0.8
        self.rib_width = 0.1
        self.x_w = x_w
        #self.F_w = F_w
        self.x_t = x_t
        self.F_t = F_t
        self.length = length
        self.t = t
        self.D = D
        self.stringers = self.make_stringers(stringer, int(n_str))
        self.Iyy = self.Iyy()
        self.max_Vz = self.max_Vz() # 2.6 * 10 **5
        self.max_My = self.max_My() # 1.4*10**6
        self.max_T = T_t
        self.J = np.pi*self.t*self.D**3/4
        self.rho = rho
        self.skin_area = self.skin_area()
        self.weight = self.get_weight()
        self.sigma_y = sigma_y
        self.buckling_factor = buckling_factor
        self.E = E
        self.poisson = poisson

    def floor_weight(self):
        floor_width = 1.9
        floor_area = floor_width * self.length
        return 5 * floor_area

    def cutout_weight(self):
        return 18 * 0.7 + 3 * 18

    def rib_weight(self):
        rib_area = np.pi*((self.D/2)**2-(self.D/2-self.rib_width)**2)
        total_area = round(self.length/self.rib_spacing) * rib_area
        return total_area*self.t*self.rho

    def skin_area(self):
        return np.pi * self.D * self.t

    def get_weight(self):
        stringer_weight = 0
        for stri in self.stringers:
            stringer_weight += stri.area * self.length * stri.rho
        bottom_skin = self.skin_area * self.length * self.rho / 2
        top_skin = self.skin_area * self.length * Glare.density / 2
        return stringer_weight + bottom_skin + top_skin + self.rib_weight() + self.floor_weight() + self.cutout_weight()

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
                stringers.append(Stringer(t, w, h, stringer.rho, stringer.sigma_y, stringer.E, stringer.poisson, angle, self.D))
            angle += str_angle
        return stringers

    def V_z(self, x):
        V = -self.weight_dist * x + self.F_w * Macaulay(x, self.x_w, 0) - self.F_t * Macaulay(x, self.x_t, 0)
        return V

    def V_z_left(self, x):
        x = self.x_w - x
        R0 = -self.weight_dist * self.x_w
        V = -R0 - self.weight_dist * x
        return V

    def V_z_right(self, x):
        R0 = -self.weight_dist*(self.length-self.x_w) - self.F_t
        V = -R0 -self.weight_dist * (x-self.x_w) - self.F_t * Macaulay(x, self.x_t, 0)
        return V

    def V_z_combined(self, x):
        if x < self.x_w:
            return self.V_z_left(x)
        else:
            return self.V_z_right(x)
    
    def max_Vz(self):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.V_z_combined(x[i]))
        return max(lst)
        
    def My(self, x):
        R0 = -(-self.weight_dist*self.length -self.F_t +self.F_w)
        M0 = self.weight_dist*self.x_w**2/2 - self.F_t * (self.x_t-self.x_w) - self.weight_dist*(self.length-self.x_w)**2/2
        print(M0)
        M = R0 * Macaulay(x, self.x_w, 1) - M0 * Macaulay(x, self.x_w, 0) - self.weight_dist * x ** 2 / 2 - self.F_t * Macaulay(x, self.x_t, 1)#+ self.F_w * Macaulay(x, self.x_w, 1)
        return M

    def My_left(self, x):
        x = self.x_w - x
        R0 = -(-self.weight_dist * (self.x_w)) #+ self.F_w
        M0 = - (self.weight_dist*(self.x_w)**2/2)
        M = R0 * Macaulay(x, 0, 1) + M0 * Macaulay(x, 0, 0) - self.weight_dist * Macaulay(x, 0, 2) / 2
        return -M

    def My_right(self, x):
        R0 = -(-self.weight_dist * (self.length -self.x_w) - self.F_t) #+ self.F_w
        M0 = (-self.F_t * (self.x_t-self.x_w) - self.weight_dist*(self.length-self.x_w)**2/2)
        M = R0 * Macaulay(x, self.x_w, 1) + M0 * Macaulay(x, self.x_w, 0) - self.weight_dist * Macaulay(x, self.x_w, 2) / 2 - self.F_t * Macaulay(x,self.x_t, 1)
        return -M

    def My_combined(self, x):
        if x < self.x_w:
            return self.My_left(x)
        else:
            return self.My_right(x)
    
    def max_My(self):
        x = np.arange(0, self.length, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.My_combined(x[i]))
        return max(lst)

    def Iyy(self):
        Iyy = 0
        for stri in self.stringers:
            Iyy += stri.Iyy
        Iyy += np.pi*self.t*self.D**3/8
        return Iyy

    def bending_stress(self, theta):
        z = np.cos(theta) * self.D/2
        return self.max_My*z/self.Iyy

    def q(self, theta):
        q = self.max_Vz/self.Iyy*np.sin(theta)*self.D**2/4*self.t
        return q

    def shear_stress(self, theta):
        return self.q(theta)/self.t + self.max_T*self.D/2/self.J

    def von_mises(self, theta):
        return np.sqrt(self.bending_stress(theta) ** 2 + 3 * self.shear_stress(theta) ** 2)

    def max_von_mises(self):
        x = np.arange(0, 2*np.pi, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.von_mises(x[i]))
        return max(lst)
    
    def max_bending_stress(self):
        x = np.arange(0, 2*np.pi, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.bending_stress(x[i]))
        return max(lst)
    
    def max_shear_stress(self):
        x = np.arange(0, 2*np.pi, 0.1)
        lst = []
        for i in range(len(x)):
            lst.append(self.shear_stress(x[i]))
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
        self.graph(self.V_z_combined, 'Shear force in z-direction', axs[0,0])
        self.graph(self.My_combined, 'Moment around the y-axis', axs[1, 0])
        self.graph_circ(self.shear_stress, 'Shear stress along the skin', axs[0, 1])
        self.graph_circ(self.bending_stress, 'Bending stress along the skin', axs[1,1])
        self.graph_circ(self.von_mises, 'Von mises stress along the skin', axs[2,1])
        #self.graph(self.momentx, 'Moment around the x-axis', axs[1, 1])
        #self.graph_displacmentz(self.graph_displacmentz, 'Displacement in z-direction', axs[2, 1])
        #self.graph(self.total_displacement_z, 'Displacement in z-direction', axs[2, 1])
        #self.graph(self.Torsiony, 'Torsional moment around the y-axis', axs[1, 2])
        #self.graph(self.total_twist, 'Twist around the y-axis', axs[2, 2])
        plt.show()

    def is_failing(self):
        return self.max_von_mises() > self.sigma_y or self.skin_buckling() or self.is_crippling()
    
    def skin_buckling(self):
        Ks = 4
        w = self.skin_width()
        tau_cr = Ks * self.E * (self.t/w)**2
        return tau_cr < self.max_shear_stress()

    def is_crippling(self):
        return self.max_bending_stress() > self.panel_crippling()
        
    def skin_width(self):
        return self.D * np.pi / len(self.stringers)
    
    def crippling(self):
        if len(self.stringers) > 0:
            return self.stringers[0].crippling(), self.max_bending_stress()
        
    def effected_width(self):
        if len(self.stringers) > 0:
            C = 4
            sigma_cc = self.stringers[0].crippling()
            return self.t/2*np.sqrt(C*np.pi**2/(12*(1-self.poisson**2)))*np.sqrt(self.E/sigma_cc)
    
    def skin_crippling(self):
        C = 4
        sigma_cc = C*np.pi**2*self.E/(12*(1-self.poisson**2))*(self.t/np.pi*self.D/2)**2
        return sigma_cc
    
    def panel_crippling(self):
        if len(self.stringers) == 0:
            return self.skin_crippling()
        effected_skin_length = 0

        A_str = 0
        sigma_cc_str = self.stringers[0].crippling()
        sigma_cc_sk = self.skin_crippling()
        for stri in self.stringers[0:len(self.stringers)//2]:
            A_str += stri.area - stri.t**2
            effected_skin_length += 2*self.effected_width()
        A_sk = (np.pi*self.D/2 - effected_skin_length)*self.t
        
        A_str += effected_skin_length*self.t
        A_sk = max(0, A_sk)
        sigma_panel = (sigma_cc_str*A_str+sigma_cc_sk*A_sk)/(A_str+A_sk)
        return sigma_panel
    
    
def create_fuselage(t_sk, n_str, material_skin, material_stringer):
    weight_ac = 84516
    n = 2.93 * 1.5
    length = wb.l_fuselage
    weight_wing = wb.weight_wing
    weight_dist = (weight_ac - weight_wing) / length * n
    t = t_sk
    D = wb.D_fuselage
    x_t = length - 0.1
    F_t = wb.L_hor * wb.b_hor * n
    x_w = wb.x_pos_wing
    F_w = weight_ac * n - weight_wing + F_t
    T_t = wb.L_ver * wb.b_vert ** 2 / 2
    buckling_factor = 0 #5 # fraction of sigma_y
    # Material
    #Stringer
    rho_str = material_stringer.density
    E_str = material_stringer.E
    sigma_y_str = material_stringer.sigma_y
    poisson_str = material_stringer.poisson
    
    #Skin
    rho_sk = material_skin.density
    E_sk = material_skin.E
    sigma_y_sk = material_skin.sigma_y
    poisson_sk = material_skin.poisson
    #G = 26.4 * 10 ** 9
    
    # Make stringer
    t_str = 0.005
    w_str = 0.04
    h_str = 0.04

    stringer = Stringer(t_str, w_str, h_str, rho_str, sigma_y_str, E_str, poisson_str)
    n_str = n_str # Has to be a multiple of 4
    
    fuselage = Fuselage(stringer, n_str, length, t, D, weight_dist, x_w, T_t, x_t, F_t, rho_sk, E_sk, sigma_y_sk, poisson_sk, buckling_factor)
    return fuselage


if __name__ == "__main__":
    t_sk = 0.003
    n_str = 20
    material_skin = AL
    material_stringer = AL7040
    fuselage = create_fuselage(t_sk, n_str, material_skin, material_stringer)
    print(fuselage.floor_weight())
    fuselage.graphs()
    print(fuselage.max_von_mises())
    print(fuselage.skin_buckling())
    print(fuselage.weight)
    print(fuselage.effected_width())
    print(fuselage.panel_crippling())