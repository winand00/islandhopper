from math import *
from ISA import script
import numpy as np
import matplotlib.pyplot as plt

C_D_0 = 0.07
C_L = 2.0

A = 7
e = 0.8
h = 2000
V = np.arange(150)
S = 72
rho = script(h)

W = 19000/2.2*9.80665

D_1 = []
for i in range(len(V)):
    D_temp = C_D_0 * 0.5 * rho * V[i]**2 * S
    D_1.append(D_temp)

D_2 = []
for i in range(len(V)):
    D_temp = ((W/(0.5*rho*V[i]**2*S))**2/(pi*A*e)) * 0.5 * rho * V[i]**2 * S
    D_2.append(D_temp)

plt.plot(V,D_1)
plt.plot(V,D_2)
plt.show()
