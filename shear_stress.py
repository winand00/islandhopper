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


def shearstress(x, z):
    if x == 0: #Left skin
        if x < h / 2:
            return -q2_z(z) + -q1_x(h/2 - z)
        else:
            return -q2_z(z) + q1_x(z-h/2)
    if z == 0: #bottom skin
        if x < w/2:
            return -q1_z(w/2-x) + -q2_x(x)
        else:
            return q3_z(w-x) + -q2_x(x)
    if abs(x-w) < 0.01: #right skins
        if x < h / 2:
            return q2_z(h-z) + -q1_x(h/2-z)
        else:
            return q2_z(h-z) + q3_x(h -z)
    if abs(z-h) < 0.01: #top skin
        if x < w / 2:
            return -q3_z(x) + q2_x(x)
        else:
            return q1_z(x-w/2) + q2_x(x)




def q1_z(s):
    return -Vz/Ixx * t * s * h/2

def q2_z(s):
    return q1_z(w/2) - Vz/Ixx * t * h/2 * s + 0.5 * s**2

def q3_z(s):
    return q2_z(h) + Vz / Ixx * t * s * h / 2


def q1_x(s):
    return Vx/Izz * t * s * w/2

def q2_x(s):
    return q1_x(h/2) + Vx/Izz * t * w/2 * s - 0.5 * s**2

def q3_x(s):
    return q2_x(w) - Vx / Izz * t * s * w / 2


N = 50

# Here are many sets of y to plot vs. x

x1 = np.arange(0, 1 + 1/N, 1/N)
y1 = np.ones(x1.shape)

x2 = np.arange(0, 1 + 1/N, 1/N)
y2 = np.zeros(x2.shape)

y3 = np.arange(0, 1 + 1/N, 1/N)
x3 = np.ones(y3.shape)

y4 = np.arange(0, 1 + 1/N, 1/N)
x4 = np.zeros(y4.shape)

x = np.hstack((x1,np.flip(x3),np.flip(x2),x4))
ys = np.hstack((y1,np.flip(y3),np.flip(y2),y4))
# We need to set the plot limits, they will not autoscale
fig, ax = plt.subplots()
ax.set_xlim(np.min(x) - 0.1, np.max(x) +0.1)
ax.set_ylim(np.min(ys)-0.1, np.max(ys)+0.1)


z = np.array([shearstress(x[i], ys[i]) for i in range(len(ys))])


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