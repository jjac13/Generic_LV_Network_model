# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a class that defines devices that combined
## can create a synthetic load profile for electric power consumed by
## households, based on their probability of being activated at a certain
## momemnt.
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0



###############################################################################
###########################   Imported modules   ##############################

import matplotlib.pyplot as plt
import seaborn as sns
from numpy import array, multiply, arange
from math import sin, acos
#import numpy as np
import csvreader

## Functions #############################################################

def Line_Currents(Node_power, Node_voltage):
        from numpy import divide
        
        I_L = []
        I_h = divide(Node_power, Node_voltage)
        n_nodes = len(Node_power)
        
        for i in range(n_nodes):
            I_L.append(0)
            for j in range(i, n_nodes):
                I_L[i] += I_h[j]
                
        return I_L
            

def Node_Voltage(V0, Z_L, I_L):
    
    V_n = []
    n_nodes = len(Z_L)
    
    for i in range(n_nodes):
        V_n.append(V0)
        for j in range(i):
            V_n[i] -= I_L[j]*Z_L[j]
    
    V_n.append(V_n[-1]-I_L[-1]*Z_L[-1])        
    return V_n


def Network_Simulation_18nodes(Node_power, V_seed, Z_L, V_error_perm = 1e-6, iteration_perm = 10):
    
    Node_power_0 = [Node_power[1], Node_power[2], (Node_power[3]+Node_power[10]), (Node_power[4]+sum(Node_power[11:15])), Node_power[5], Node_power[6]+Node_power[15], Node_power[7],  Node_power[8], Node_power[9]+Node_power[16],  Node_power[10]+Node_power[17]]
#    V_error = 200
    iteration = 0
    
    while iteration < iteration_perm:
        # First iteration
        V_seed_0 = V_seed[0:10]
        I_L_0 = Line_Currents(Node_power_0, V_seed_0)
        V_n_0 = Node_Voltage(V_seed[0], Z_L[0:10], I_L_0)
        
        ## Analyzing the branches
        ## Branch 3-11
        #    Node_power_3_11 = Node_power[10]
        I_L_3_11 = I_L_0[2] - I_L_0[3]
        V_n_3_11 = V_n_0[3] - I_L_3_11 * Z_L[10]
        
        ## Branch 4-15
        Node_power_4_15 = Node_power[11:15]
        I_L_4_15 = Line_Currents(Node_power_4_15, V_n_0[4])
        V_n_4_15 = Node_Voltage(V_n_0[4], Z_L[11:15], I_L_4_15)
        
        ## Branch 6-16
        #    Node_power_6_16 = Node_power[15]
        I_L_6_16 = I_L_0[5] - I_L_0[6]
        V_n_6_16 = V_n_0[6] - I_L_6_16 * Z_L[15]
        
        ## Branch 9-17
        #    Node_power_9_17 = Node_power[16]
        I_L_9_17 = I_L_0[8] - I_L_0[9]
        V_n_9_17 = V_n_0[6] - I_L_9_17 * Z_L[16]
        
        ## Branch 10-18
        #    Node_power_10_18 = Node_power[17]
        I_L_10_18 = I_L_0[9]
        V_n_10_18 = V_n_0[9] - I_L_10_18 * Z_L[17]
        
        
        V_n = V_n_0 + [V_n_3_11] + V_n_4_15[1:5] + [V_n_6_16] + [V_n_9_17] + [V_n_10_18]
        I_L = I_L_0 + [I_L_3_11] + I_L_4_15[1:5] + [I_L_6_16] + [I_L_9_17] + [I_L_10_18]
        
        V_error = PolarComponents(sum([PolarComponents(i - j) for i,j in zip(V_seed, V_n)]))
        V_seed = V_n
        iteration +=1
#        print('The error of the iteration: ', iteration, ' is: ', V_error)
        
        if abs(V_error) < V_error_perm:
            break
    
    return [V_n, I_L]

def Network_Simulation_6nodes(Node_power, V_seed, Z_L, V_error_perm = 1e-6, iteration_perm = 10):
    
    Node_power_0 = [Node_power[1], Node_power[2], Node_power[3], Node_power[4], Node_power[5], Node_power[5]]
#    V_error = 200
    iteration = 0
    
    while iteration < iteration_perm:
        # First iteration
        V_seed_0 = V_seed[0:6]
        I_L_0 = Line_Currents(Node_power_0, V_seed_0)
        V_n_0 = Node_Voltage(V_seed[0], Z_L[0:6], I_L_0)
                    
        V_n = V_n_0
        I_L = I_L_0
        
        V_error = PolarComponents(sum([PolarComponents(i - j) for i,j in zip(V_seed, V_n)]))
        V_seed = V_n
        iteration +=1
#        print('The error of the iteration: ', iteration, ' is: ', V_error)
        
        if abs(V_error) < V_error_perm:
            break
    
    return [V_n, I_L]

def PolarComponents(complex_number, magnitude = True, rad = False):
    from cmath import polar, pi
    
    [modulus, angle] = polar(complex_number)
    
    if magnitude == True:
#        return modulus
        if complex_number.real > 0:
            return modulus
        else:
            return -modulus
    else:
        if rad == False:
            return 180*angle/pi    
        else:
            return angle
                
    
        
## Main code #############################################################

V0 = 400                    # Voltage at the bulk power source.
n_nodes = 18                # Number of nodes.
V_seed = [V0]*n_nodes      # Seed voltage for the nodes.
#Z_unit = (0.22)/1000  # Impedance per meter of distribution wire.
#Z_unit = (0.22+0.37j)/1000  # Impedance per meter of distribution wire.
#Z_unit = [0, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 0.49+0.471j, 2.334 + 1.545j, 0.733 + 0.57j, 0.733 + 0.57j, 0.733 + 0.57j, 0.733 + 0.57j, 1.285 + 0.865j, 2.334 + 1.454j, 1.926 + 1.265j]
# 18-node
Z_unit = [0, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 0.162 +0.07j, 1.539+0.076j, 0.265+0.07065j, 0.265+0.07065j, 0.265+0.07065j, 0.265+0.07065j, 0.229+.0719j, 1.539+0.076j, 1.113+0.0735j]
distances_z = [35, 35, 35, 35, 35, 35, 35, 35, 35,  35, 30, 35, 35, 35, 30, 30, 30, 30]  # Distances between nodes.

# 6-node
#Z_unit = [0, 0.549 +0.072j, 0.549 +0.072j, 0.549 +0.072j, 0.549 +0.072j, 0.549 +0.072j]
#distances_z = [100]*6  # Distances between nodes.

Z_L = multiply(Z_unit, distances_z)
Z_L = Z_L/1000
    
PowerData = csvreader.read_data(csv='Case_3_Summer_PVT.csv', address='')
PowerData.data2rows()


### 18-node circuit    
#Node_power = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.57 + 3.45j, 0, 0, 0, (17.04-0) + 10.56j, (17.4-0) + 10.56j, 5.57 + 3.45j,  (8.2-0) + 5.08j]
#Node_power = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.57 + 3.45j, 0, 0, 0, (17.04-5.5) + 10.56j, (17.4-4) + 10.56j, 5.57 + 3.45j,  (8.2-3) + 5.08j]
Node_power_c1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.93 + 6.15j, 0, 0, 0, 9.93 + 6.15j, 9.93 + 6.15j, 9.93 + 6.15j,  9.93 + 6.15j]
Node_power_c1_HP = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (9.93+5.4) + 6.15j, 0, 0, 0, (9.93+8.1) + 6.15j, (9.93+16.2) + 6.15j, (9.93+16.2) + 6.15j,  (9.93+16.2) + 6.15j]
#Node_power = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 42.5, 0, 0, 0, 42.5, 42.5, 42.5,  42.5]
Node_power_c1 = [i*1000 for i in Node_power_c1]
Node_power_c1_HP = [i*1000 for i in Node_power_c1_HP]

## 6-node circuit    
#Node_power_c1 = [0, 0, 0, 0, 0, 0]
#Node_power_c1_HP = [0, 6, 6, 6, 6, 6]
#Node_power_c1 = [i*1000 for i in Node_power_c1]
#Node_power_c1_HP = [i*1000 for i in Node_power_c1_HP]

Total_currents = []
Total_voltages = []

Current_mag = []
Current_phase = []
Voltage_mag = []
Voltage_phase = []

#Total_power = [sum(power) for power in PowerData.row]


[V_n_c1, I_L_c1] = Network_Simulation_18nodes(Node_power_c1, V_seed, Z_L)
[V_n_c1_HP, I_L_c1_HP] = Network_Simulation_18nodes(Node_power_c1_HP, V_seed, Z_L)

#for Power_timestep in PowerData.row:
##    Node_power = [(i + (i*sin(acos(0.85))/0.85)*1j)*1000 for i in Power_timestep]   # To watts, and to 6 houses
#    Node_power = [i*1000 for i in Power_timestep]   # To watts, and to 6 houses    
#
#    I_L = Line_Currents(Node_power, V_seed)
#    V_n = Node_Voltage(V0, Z_L, I_L)
#    
#    Total_currents.append(I_L)
#    Current_mag.append([PolarComponents(current) for current in I_L])
#    Current_phase.append([PolarComponents(current, False) for current in I_L])
#    
#    Total_voltages.append(V_n)
#    Voltage_mag.append([PolarComponents(voltage)/PolarComponents(V_n[0]) for voltage in V_n])
#    Voltage_phase.append([PolarComponents(voltage, False) for voltage in V_n])
#    
#    V_seed = V_n[1:n_nodes+1]


Current_mag_c1 = [PolarComponents(current) for current in I_L_c1]
Current_phase_c1 = [PolarComponents(current, False) for current in I_L_c1]
Voltage_mag_c1 = [PolarComponents(voltage)/PolarComponents(V_n_c1[0]) for voltage in V_n_c1]
Voltage_phase_c1 = [PolarComponents(voltage, False) for voltage in V_n_c1]


Current_mag_c1_HP = [PolarComponents(current) for current in I_L_c1_HP]
Current_phase_c1_HP = [PolarComponents(current, False) for current in I_L_c1_HP]
Voltage_mag_c1_HP = [PolarComponents(voltage)/PolarComponents(V_n_c1[0]) for voltage in V_n_c1_HP]
Voltage_phase_c1_HP = [PolarComponents(voltage, False) for voltage in V_n_c1_HP]



Letter_size = 16

plt.rc('font', size=Letter_size)          # controls default text sizes
plt.rc('axes', titlesize=Letter_size)     # fontsize of the axes title
plt.rc('axes', labelsize=Letter_size)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=Letter_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=Letter_size)    # fontsize of the tick labels
plt.rc('legend', fontsize=Letter_size)    # legend fontsize
plt.rc('figure', titlesize=Letter_size)  # fontsize of the figure title

#fig1 = plt.figure(1)
#plt.plot(Total_power, 'k', label='Total power demanded in the circuit')
#plt.xlim([0, len(Total_power)])
#plt.xlabel('Timestep')
#plt.ylabel('Power [kW]')
#plt.title('Total power demanded in the circuit')
##plt.legend(loc='upper right')
#plt.grid()
#plt.show()

#fig2 = plt.figure(2)
#plt.plot(Current_mag_c1, 'b', label='Line current out of the nodes during summer')
#plt.plot(Current_mag_c1_HP, 'r', label='Line current out of the nodes during winter')
#plt.xlim([0, 10])
#plt.xlabel('Nodes')
#plt.ylabel('Current [A]')
#plt.title('Line current out of the nodes')
#plt.legend(loc='upper right')
#plt.grid()
#plt.show()
#
fig3 = plt.figure(3)
#plt.plot(Voltage_mag_c1, 'b', label='Voltage in the node during summer')
plt.plot(Voltage_mag_c1_HP, 'r')#, label='Voltage in the node during winter')
plt.xlim([1, 6])
plt.xticks(arange(1, 7, step=1))  # Set label locations. 
plt.xlabel('Nodes')
#plt.ylim([0.965, 1.001])
plt.ylabel('Voltage [pu]') # 'Voltage [V]'
#plt.title('Voltage per node')
#plt.legend(loc='lower left')
plt.grid()
plt.show()

#fig4 = plt.figure(4)
#sns.heatmap(Voltage_mag, cmap='rainbow', cbar_kws={'label': 'Voltage [pu]'})
#plt.xlabel('Nodes')
#plt.ylabel('Hour of the day')
##plt.title('Voltage in the node')
##plt.legend(loc='upper right')
##plt.grid()
#plt.show()


#fig5 = plt.figure(5)
#sns.heatmap(Current_mag, cmap='rainbow', cbar_kws={'label': 'Current [A]'},) # annot=True, fmt=".01f", 
#plt.xlabel('Nodes')
#plt.ylabel('Hour of the day')
#plt.title('Line current out of the nodes')
##plt.legend(loc='upper right')
##plt.grid()
#plt.show()

#fig6 = plt.figure(6)
#plt.plot(Voltage_mag[32], 'g', label='Voltage in the node at 8:00 am')
#plt.plot(Voltage_mag[48], 'r', label='Voltage in the node at 12:00 am')
#plt.plot(Voltage_mag[84], 'b', label='Voltage in the node at 09:00 pm')
#plt.xlim([1, 18])
#plt.xticks(arange(1, 19, step=1))  # Set label locations. 
#plt.xlabel('Node')
#plt.ylabel('Voltage [pu]')
##plt.title('Voltage in the node')
#plt.legend(loc='upper left')
#plt.grid()
#plt.show()