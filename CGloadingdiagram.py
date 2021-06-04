import numpy as np
import matplotlib.pyplot as plt
import math as m
from matplotlib.colors import *
from Class2_weight_estimation import cg_OEW,W_OEW

'''
Function 'loading_diagram' can be used to plot the loading diagram of either the RJ85 or the RJXX.
Inputs of function are:
: Hopper (input as a string)f
plot: single or combined (input as a string)
'''
def cgboundaries():
    min_xcg_Hopper = min(loading_diagram()[0]) - 0.02
    max_xcg_Hopper = max(loading_diagram()[2]) + 0.02
    return min_xcg_Hopper, max_xcg_Hopper

def loading_diagram():
    # Define Hopper input parameters
    MAC                     = 2    #??                                          # [m]
    MAC_start               = 7   # ??                                          # [m]
    OEW                     = W_OEW    # ??                                          # from class2 script [kg]
    X_oew                   = cg_OEW                                           # [percentage of MAC] (Still Assumed 0.3569)!
    M_fuel                  = 110                                              # To not exceed the MTOW, a fuel weight of 6249 [kg] is used!!!                                              # [-]
    M_pax                   = 82                                               # [kg] 77kg of passenger + 5kg carry on baggage
    seat_pitch              = 0.762     #30 inches (ADSEE)
    nosecone_length          = 3.056                         # ends at middle of row 1
    starting_length         = nosecone_length - 0.762     # 1 seat distance before middle of row 1
    X_fuel                  = 6.5         # ??                                          # [m] location of cg fuel                                            # [m] location of cg tank
    M_front_cargo           = 95                     # [kg]
    M_rear_cargo            = 95                      # [kg]
    x_front_cargo_absolute  = 3    # ??
    x_rear_cargo_absolute   = 8    # ??
    X_front_cargo           = (x_front_cargo_absolute - MAC_start)/MAC
    X_rear_cargo            = (x_rear_cargo_absolute - MAC_start)/MAC
    emergency_space         = 0.15       # [m]
    # ____________________________________________________________________________

    # Function to calculate the cg by adding extra mass:
    def cg_calc(CurrentMass,Current_cg,Extra_mass,cg_extra_mass):
        New_mass = CurrentMass + Extra_mass
        New_cg = (Current_cg*CurrentMass+Extra_mass*cg_extra_mass)/New_mass
        return New_mass, New_cg

    cg_pos_front_load  = [X_oew]
    weight_front_load = [OEW]


    # First load the cargo from the front :
    New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],M_front_cargo,X_front_cargo)
    cg_pos_front_load.append(New_cg)
    weight_front_load.append(New_mass)
    # Now add the cargo to the back:
    New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],M_rear_cargo,X_rear_cargo)
    cg_pos_front_load.append(New_cg)
    weight_front_load.append(New_mass)


    #First load the passengers from the front, 2 window rows first
    for i in range(8):
        if i>0 and i<4: # First 3 rows
            New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],2*M_pax,((starting_length+i*seat_pitch-MAC_start)/MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)
        if i>3 and i<7: # Next 3 rows
            New_mass, New_cg = cg_calc(weight_front_load[-1], cg_pos_front_load[-1], 2 * M_pax,
                                       ((starting_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)
        if i>6: #Last row 1 person
            New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],1*M_pax,((starting_length+i*seat_pitch+emergency_space-MAC_start)/MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)

    #Load the passengers from the front, 1 middle row seaters now
    for i in range(7):
        if i>0 and i<4: #First 3 rows
            New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],1*M_pax,((starting_length+i*seat_pitch-MAC_start)/MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)
        if i > 3:  # Next 3 rows
            New_mass, New_cg = cg_calc(weight_front_load[-1], cg_pos_front_load[-1], 1 * M_pax,
                                       ((starting_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)
    # ____________________________________________________________________________

    cg_pos_rear_load  = [X_oew]
    weight_rear_load = [OEW]

    # This time round, we load everything from the rear of the plane, starting with cargo:
    New_mass, New_cg = cg_calc(weight_rear_load[-1],cg_pos_rear_load[-1],M_rear_cargo,X_rear_cargo)
    cg_pos_rear_load.append(New_cg)
    weight_rear_load.append(New_mass)
    New_mass, New_cg = cg_calc(weight_rear_load[-1],cg_pos_rear_load[-1],M_front_cargo,X_front_cargo)
    cg_pos_rear_load.append(New_cg)
    weight_rear_load.append(New_mass)

    # Load the passengers from the rear, Last row only 1 person extra
    for i in range(8):
        if i > 6:  # Only last row
            New_mass, New_cg = cg_calc(weight_rear_load[-1], cg_pos_rear_load[-1], 1 * M_pax,
                                       ((starting_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_rear_load.append(New_cg)
            weight_rear_load.append(New_mass)

    # Load window passengers from rear
    for i in range(6, 0, -1):
        if i > 3 and i < 7: # Row 4-6
            New_mass, New_cg = cg_calc(weight_rear_load[-1], cg_pos_rear_load[-1], 2 * M_pax,
                                       ((starting_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_rear_load.append(New_cg)
            weight_rear_load.append(New_mass)
        if i < 4: # Row 1-3
            New_mass, New_cg = cg_calc(weight_rear_load[-1], cg_pos_rear_load[-1], 2 * M_pax,
                                       ((starting_length + i * seat_pitch - MAC_start) / MAC))
            cg_pos_rear_load.append(New_cg)
            weight_rear_load.append(New_mass)

    # Load middle passengers from rear
    for i in range(6, 0, -1):
        if i > 3 and i < 7:  # Row 4-6
            New_mass, New_cg = cg_calc(weight_rear_load[-1], cg_pos_rear_load[-1], 1 * M_pax,
                                       ((starting_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_rear_load.append(New_cg)
            weight_rear_load.append(New_mass)
        if i < 4:  # Row 1-3
            New_mass, New_cg = cg_calc(weight_rear_load[-1], cg_pos_rear_load[-1], 1 * M_pax,
                                       ((starting_length + i * seat_pitch - MAC_start) / MAC))
            cg_pos_rear_load.append(New_cg)
            weight_rear_load.append(New_mass)

    # ____________________________________________________________________________

    # Tanking the fuel, since there is only one way to do this it is both from the front and rear
    New_mass, New_cg = cg_calc(weight_rear_load[-1],cg_pos_rear_load[-1],M_fuel, (X_fuel- MAC_start) / MAC)
    cg_pos_rear_load.append(New_cg)
    weight_rear_load.append(New_mass)


    return [cg_pos_front_load, weight_front_load, cg_pos_rear_load, weight_rear_load]

min_xcg_Hopper = min(loading_diagram()[0]) - 0.02
max_xcg_Hopper = max(loading_diagram()[2]) + 0.02

def plot_loadings():
    fig, ax = plt.subplots()
    plt.plot(loading_diagram()[0][0:3]  , loading_diagram()[1][0:3]  , label = 'Cargo', color='b')
    plt.plot(loading_diagram()[0][2:10] , loading_diagram()[1][2:10] , label = 'Window seats', color='r')
    plt.plot(loading_diagram()[0][9:16], loading_diagram()[1][9:16], label = 'Middle seats', color='c')
    plt.plot(loading_diagram()[2][-2:]  , loading_diagram()[3][-2:]  , label = 'Fuel', color='orange')
    plt.plot(loading_diagram()[2][0:3]  , loading_diagram()[3][0:3]  , color='b')
    plt.plot(loading_diagram()[2][2:10] , loading_diagram()[3][2:10] , color='r')
    plt.plot(loading_diagram()[2][9:16], loading_diagram()[3][9:16],  color='c')

    plt.title('Loading Diagram Hopper')
    plt.xlabel(r'$X_{cg}/MAC$ [-]')
    plt.ylabel('Weight [kg]')
    plt.ylim((5000,8700))
    plt.xlim(-0.05,0.5)
    plt.xticks(np.arange(0.05,0.50,0.05))
    ax.set(facecolor='w')
    plt.axvline(min_xcg_Hopper, ymin=0, ymax=1, color='black', linestyle='-')
    plt.axvline(max_xcg_Hopper, ymin=0, ymax=1, color='black', linestyle='-')
    # Get 2% margin lines:

    # Legend:
    plt.legend(loc=2, prop={'size': 8})
    plt.grid(axis='y')
    plt.show()

if __name__ == "__main__":
    print('Minimum x_cg Hopper:', min_xcg_Hopper)
    print('Maximum x_cg Hopper:', max_xcg_Hopper)
    print()

    plot_loadings()
