import numpy as np
import matplotlib.pyplot as plt
import math as m
from matplotlib.colors import *
'''
Function 'loading_diagram' can be used to plot the loading diagram of either the RJ85 or the RJXX.
Inputs of function are:
aircraft: RJ85 or RJXX (input as a string)
plot: single or combined (input as a string)
'''


def loading_diagram(aircraft):
    # Define RJ85 input parameters
    if aircraft=='Hopper':
        MTOW                    = 8618.25503                                             # [kg]
        MAC                     = 2    #??                                          # [m]
        MAC_start               = 7   # ??                                          # [m]
        OEW                     = 0.6*MTOW    # ??                                          # [kg]
        X_oew                   = 0.3569                                           # [percentage of MAC] Still Assumed!
        M_fuel                  = 6350                                              # To not exceed the MTOW, a fuel weight of 6249 [kg] is used!!!
        M_payload               = 1748   #??                                            # [kg] both cargo and passengers
        N_pax                   = 19                                               # [-]
        M_pax                   = 77                                               # [kg]
        N_seat_abreast          = 3                                                 # [-]
        N_rows                  = 6
        seat_pitch              = 0.762     #??
        cabin_length            = 7.366  #Jetstream                                           # [m]
        length                  = 16.38  #Jetstream + 2 [m]
        tailcone_length         = 4.2084                                             #3/5*(length-2m-cabin_length)
        nosecone_length         = 2.8056                                             #2/5*(length-2m-cabin_length)
        cabin_start             = nosecone_length
        cabin_start_relative    = (nosecone_length-MAC_start)/MAC
        cabin_end_relative      = (nosecone_length+cabin_length-MAC_start)/MAC
        effective_cabin_length  = cabin_length-7*seat_pitch
        galley_length           = (cabin_length-effective_cabin_length)/2
        X_fuel_absolute         = MAC/2                                             # Schatting, nog niet uitgerekend!!!!!!!!
        X_fuel                  = 0.6                              # [percentage of the MAC]
        X_tank                  = 0.4                                               # [percentage of the MAC]
        M_front_cargo           = 285    # No cargo spaces                  # [kg]
        M_rear_cargo            = 0                      # [kg]
        x_front_cargo_absolute  = 2    # ??
        x_rear_cargo_absolute   = 0
        X_front_cargo           = (x_front_cargo_absolute - MAC_start)/MAC 
        X_rear_cargo            = (x_rear_cargo_absolute - MAC_start)/MAC
        emergency_space         = 0.2    # ??   # [m]
    # ____________________________________________________________________________

    # Define RJXX input parameters
    '''
    OEW decreases by 10%, MTOW remains the same. The weight reduction in OEW will be compensated by
    an equal increase in MPW (max. payload), by means of extra cargo (no extra passengers). 75% of 
    the extra cargo will be stored in the front cargo hold, 25% will be stored in the rear cargo hold.
    In addition, X_oew shifts 60 [cm] (0.6 [m]) forward. 
    '''
    if aircraft=='RJXX':
        MTOW                    = 42184                                             # [kg]
        MAC                     = 3.17                                              # [m] (given from aircraft data)
        MAC_start               = 10.60                                             # [m]
        OEW                     = 0.90 * 24820                                      # [kg]
        X_oew_RJ85              = 0.3569                                           # [% MAC] (still assumed!!!)
        X_oew_absolute          = X_oew_RJ85 * MAC + MAC_start - 0.6                # [m]
        X_oew                   = 0.1703                                            # [% MAC]
        M_fuel                  = 6350                                              # [kg] 
        X_fuel_absolute         = MAC/2                                             # Schatting, nog niet uitgerekend!!!!!!!!
        X_fuel                  = 0.6                               # [percentage of the MAC]
        X_tank                  = 0.4                                               # [percentage of the MAC]
        M_payload               = 11014 + 0.1 * 24820                               # [kg] both cargo and passengers
        N_pax                   = 100                                               # [-]
        M_pax                   = 95                                             # [kg]
        N_seat_abreast          = 6                                                 # [-]
        N_rows                  = N_pax/N_seat_abreast
        seat_pitch              = 0.762
        cabin_length            = 17.11                                             # [m]
        length                  = 28.55
        tailcone_length         = 7.123                                             #3/5*(length-cabin_length)
        nosecone_length         = 5.7                                             #2/5*(length-cabin_length)
        cabin_start             = nosecone_length
        cabin_start_relative    = (nosecone_length-MAC_start)/MAC
        cabin_end_relative      = (nosecone_length+cabin_length-MAC_start)/MAC
        effective_cabin_length  = cabin_length-17*seat_pitch
        galley_length           = (cabin_length-effective_cabin_length)/2

        # Add extra cargo weight (75% front, 25% rear):
        M_front_cargo           = 757  + 0.75*0.10*24820        # [kg]
        M_rear_cargo            = 757  + 0.25*0.10*24820         # [kg]
        x_front_cargo_absolute  = 6.8
        x_rear_cargo_absolute   = 17.92
        X_front_cargo           = (x_front_cargo_absolute - MAC_start)/MAC 
        X_rear_cargo            = (x_rear_cargo_absolute - MAC_start)/MAC
        emergency_space = 0.2  # ??   # [m]
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
    for i in range(4):
        print(i)
        if i>0: #Otherwise it also takes into account an extra 2 passengers when i=0
            New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],2*M_pax,((nosecone_length+i*seat_pitch-MAC_start)/MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)

    #Load the passengers from the front, 1 middle row seaters now
    for i in range(4):
        if i>0: #Otherwise it also takes into account an extra 2 passengers when i=0
            New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],1*M_pax,((nosecone_length+i*seat_pitch-MAC_start)/MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)

    # Load the passengers from the front, 1 middle row seaters now
    for i in range(7):
        if i > 3:  # Otherwise it also takes into account an extra 2 passengers when i=0
            New_mass, New_cg = cg_calc(weight_front_load[-1], cg_pos_front_load[-1], 1 * M_pax,
                                       ((nosecone_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)

    #Load the passengers from the front, Last row only 1 person extra
    for i in range(8):
        if i>6: #Otherwise it also takes into account an extra passenger when i=0
            New_mass, New_cg = cg_calc(weight_front_load[-1],cg_pos_front_load[-1],1*M_pax,((nosecone_length+i*seat_pitch+emergency_space-MAC_start)/MAC))
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
        if i > 6:  # Otherwise it also takes into account an extra passenger when i=0
            New_mass, New_cg = cg_calc(weight_front_load[-1], cg_pos_front_load[-1], 1 * M_pax,
                                       ((nosecone_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
            cg_pos_front_load.append(New_cg)
            weight_front_load.append(New_mass)

    # Load passengers from rear, row 4-6
    for i in range(6, 3, -1):
        New_mass, New_cg = cg_calc(weight_front_load[-1], cg_pos_front_load[-1], 1 * M_pax,
                                   ((nosecone_length + i * seat_pitch + emergency_space - MAC_start) / MAC))
        cg_pos_front_load.append(New_cg)
        weight_front_load.append(New_mass)

    # Load the passengers from the rear, row 1-3
    for i in range(3, 0, -1):
        New_mass, New_cg = cg_calc(weight_front_load[-1], cg_pos_front_load[-1], 1 * M_pax,
                                   ((nosecone_length + i * seat_pitch - MAC_start) / MAC))
        cg_pos_front_load.append(New_cg)
        weight_front_load.append(New_mass)
    # ____________________________________________________________________________

    # Tanking the fuel, since there is only one way to do this it is both from the front and rear
    New_mass, New_cg = cg_calc(weight_rear_load[-1],cg_pos_rear_load[-1],M_fuel,X_fuel)
    cg_pos_rear_load.append(New_cg)
    weight_rear_load.append(New_mass)

    # New_mass, New_cg = cg_calc(weight_rear_load[-1],cg_pos_rear_load[-1],M_fuel,X_fuel)
    # cg_pos_rear_load.append(New_cg)
    # weight_rear_load.append(New_mass)

    # plt.plot(cg_pos_front_load,weight_front_load, label = 'Forward Loading')
    # plt.plot(cg_pos_rear_load,weight_rear_load, label = 'Rearward Loading')

    return [cg_pos_front_load, weight_front_load, cg_pos_rear_load, weight_rear_load]


min_xcg_RJ85 = min(loading_diagram('Hopper')[0]) - 0.02
min_xcg_RJXX = min(loading_diagram('RJXX')[0]) - 0.02
max_xcg_RJ85 = max(loading_diagram('Hopper')[2]) + 0.02
max_xcg_RJXX = max(loading_diagram('RJXX')[2]) + 0.02

'''
Function inputs:
- 'RJ85'
- 'RJXX'
- 'Both'
'''
def plot_loadings(aircraft):
    MTOW = 42184

    if aircraft=='Hopper':
        fig, ax = plt.subplots()
        plt.plot(loading_diagram(aircraft)[0][0:3]  , loading_diagram(aircraft)[1][0:3]  , label = 'Cargo', color='b')
        plt.plot(loading_diagram(aircraft)[0][2:20] , loading_diagram(aircraft)[1][2:20] , label = 'Window seats', color='r')
        plt.plot(loading_diagram(aircraft)[0][19:37], loading_diagram(aircraft)[1][19:37], label = 'Middle seats', color='c')
        plt.plot(loading_diagram(aircraft)[0][36:54], loading_diagram(aircraft)[1][36:54], label = 'Aisle seats', color='k')
        plt.plot(loading_diagram(aircraft)[2][-2:]  , loading_diagram(aircraft)[3][-2:]  , label = 'Fuel', color='orange')
        plt.plot(loading_diagram(aircraft)[2][0:3]  , loading_diagram(aircraft)[3][0:3]  , color='b')
        plt.plot(loading_diagram(aircraft)[2][2:20] , loading_diagram(aircraft)[3][2:20] , color='r')
        plt.plot(loading_diagram(aircraft)[2][19:37], loading_diagram(aircraft)[3][19:37],  color='c')
        plt.plot(loading_diagram(aircraft)[2][36:54], loading_diagram(aircraft)[3][36:54],  color='k')
        
        plt.title('Loading Diagram Bae AVRO RJ-85')
        plt.xlabel(r'$X_{cg}/MAC$ [-]')
        plt.ylabel('Weight [kg]')
        plt.ylim((5000,9000))
        plt.xlim(0.05,0.9)
        plt.xticks(np.arange(0.20,0.70,0.05))
        ax.set(facecolor='w')
        plt.axvline(min_xcg_RJ85, ymin=0, ymax=1, color='black', linestyle='-')
        plt.axvline(max_xcg_RJ85, ymin=0, ymax=1, color='black', linestyle='-')
        # Get 2% margin lines:
        
        # Legend:
        plt.legend(loc=2, prop={'size': 8})
        plt.grid(axis='y')
        plt.show()


    if aircraft=='RJXX':
        fig, ax = plt.subplots()
        plt.plot(loading_diagram(aircraft)[0][0:3]  , loading_diagram(aircraft)[1][0:3]  , label = 'Cargo', color='b')
        plt.plot(loading_diagram(aircraft)[0][2:20] , loading_diagram(aircraft)[1][2:20] , label = 'Window seats', color='r')
        plt.plot(loading_diagram(aircraft)[0][19:37], loading_diagram(aircraft)[1][19:37], label = 'Middle seats', color='c')
        plt.plot(loading_diagram(aircraft)[0][36:54], loading_diagram(aircraft)[1][36:54], label = 'Aisle seats', color='k')
        plt.plot(loading_diagram(aircraft)[2][-2:]  , loading_diagram(aircraft)[3][-2:]  , label = 'Fuel ', color='orange')
        plt.plot(loading_diagram(aircraft)[2][0:3]  , loading_diagram(aircraft)[3][0:3]  ,  color='b')
        plt.plot(loading_diagram(aircraft)[2][2:20] , loading_diagram(aircraft)[3][2:20] ,  color='r')
        plt.plot(loading_diagram(aircraft)[2][19:37], loading_diagram(aircraft)[3][19:37],  color='c')
        plt.plot(loading_diagram(aircraft)[2][36:54], loading_diagram(aircraft)[3][36:54],  color='k')

        plt.title('Loading Diagram Bae AVRO RJ-XX')
        plt.xlabel(r'$X_{cg}/MAC$ [-]')
        plt.ylabel('Weight [kg]')
        plt.ylim((20000,45000))
        plt.xlim(-0.1,0.4)
        plt.xticks(np.arange(-0.1,0.4,0.05))

        ax.set(facecolor='w')
        plt.axvline(min_xcg_RJXX, ymin=0, ymax=1, color='black', linestyle='-')
        plt.axvline(max_xcg_RJXX, ymin=0, ymax=1, color='black', linestyle='-')
        # Get 2% margin lines:
       
        # Legend:
        plt.legend(loc=2, prop={'size': 8})
        plt.grid(axis='y')
        plt.show()


    if aircraft=='Both':
        fig, ax = plt.subplots()
        plt.plot(loading_diagram('RJ85')[0][0:3]  , loading_diagram('RJ85')[1][0:3]   , color='b', linestyle='--')
        plt.plot(loading_diagram('RJ85')[0][2:20] , loading_diagram('RJ85')[1][2:20]  , color='r', linestyle='--')
        plt.plot(loading_diagram('RJ85')[0][19:37], loading_diagram('RJ85')[1][19:37] , color='c', linestyle='--')
        plt.plot(loading_diagram('RJ85')[0][36:54], loading_diagram('RJ85')[1][36:54] , color='k', linestyle='--')
        plt.plot(loading_diagram('RJ85')[2][-2:]  , loading_diagram('RJ85')[3][-2:]   , color='orange', linestyle='--')
        plt.plot(loading_diagram('RJ85')[2][0:3]  , loading_diagram('RJ85')[3][0:3]   , color='b', linestyle='--')
        plt.plot(loading_diagram('RJ85')[2][2:20] , loading_diagram('RJ85')[3][2:20]  , color='r', linestyle='--')
        plt.plot(loading_diagram('RJ85')[2][19:37], loading_diagram('RJ85')[3][19:37] , color='c', linestyle='--')
        plt.plot(loading_diagram('RJ85')[2][36:54], loading_diagram('RJ85')[3][36:54] , color='k', linestyle='--')

        plt.plot(loading_diagram('RJXX')[0][0:3]  , loading_diagram('RJXX')[1][0:3]  , color='b')
        plt.plot(loading_diagram('RJXX')[0][2:20] , loading_diagram('RJXX')[1][2:20] , color='r')
        plt.plot(loading_diagram('RJXX')[0][19:37], loading_diagram('RJXX')[1][19:37], color='c')
        plt.plot(loading_diagram('RJXX')[0][36:54], loading_diagram('RJXX')[1][36:54], color='k')
        plt.plot(loading_diagram('RJXX')[2][-2:]  , loading_diagram('RJXX')[3][-2:]  , color='orange')
        plt.plot(loading_diagram('RJXX')[2][0:3]  , loading_diagram('RJXX')[3][0:3]  , color='b')
        plt.plot(loading_diagram('RJXX')[2][2:20] , loading_diagram('RJXX')[3][2:20] , color='r')
        plt.plot(loading_diagram('RJXX')[2][19:37], loading_diagram('RJXX')[3][19:37], color='c')
        plt.plot(loading_diagram('RJXX')[2][36:54], loading_diagram('RJXX')[3][36:54], color='k')

        plt.title('Loading Diagram of the Bae AVRO RJ-85 and Bae Avro RJ-XX')
        plt.xlabel(r'$X_{cg}/MAC$ [-]')
        plt.ylabel('Weight [kg]')
        plt.ylim((21000,45000))
        plt.xticks(np.arange(-0.05,0.50,0.05))
        plt.grid(axis='y')
        ax.set(facecolor='w')
        
        lines = ax.get_lines()
        legend1 = plt.legend([lines[i] for i in [9,10,11,12,13]], ['Cargo', 'Window seats', 'Middle seats', 'Aisle seats', 'Fuel'], loc=2)
        ax.add_artist(legend1)

        plt.show()  

print('Minimum x_cg RJ85:', min_xcg_RJ85)
print('Maximum x_cg RJ85:', max_xcg_RJ85)
print()
print('Minimum x_cg RJXX:', min_xcg_RJXX)
print('Maximum x_cg RJXX:', max_xcg_RJXX)



plot_loadings('Hopper')
plot_loadings('RJXX')
plot_loadings('Both')