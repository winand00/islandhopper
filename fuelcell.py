import numpy as np
import matplotlib.pyplot as plt
#inputs

#external
P_cruise = 500000 #Watt
P_take_off = 1300000  #watt

#constants
# y = ax+b : y = gewicht, x = aantal cellen

fuel_cell_a = 0.0725
fuelcell_b = 9.37
efficiency_constant = 1.482

A = np.linspace(0,450,45)

def voltage(A):
    return -0.047*np.log(A) + 0.9782

def efficiency(V):
    return V/efficiency_constant

V = voltage(A)
P = V*A

plt.plot(A, V)
plt.show()
plt.plot(A,P)
plt.show()
plt.plot(A, efficiency(V))
plt.show()