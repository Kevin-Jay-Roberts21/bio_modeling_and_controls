# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:12:18 2024

@author: Kevin Roberts
"""

# Note: the results are output to the file "output.txt"

from matplotlib import pyplot as plt
import math
import numpy as np

# Non-Dimensional Constants
alpha = 0.01
beta = 0.1
Dc = 25
mu = 13786
delta = 110
A = 560
sigma = 4000
c_hat = 3

# dt = .5 * (dx*2) / Dc
dt = 1/12
dx_max = math.sqrt(2 * Dc * dt) # 2.0412414523193148
dx = .5

# Supplementary Equations
# n = cell density (nondimensional) and c = EGF concentration (nondimensional)
def B(n):                   
    if n <= .2:
        return sigma
    elif n >= .4:
        return 0
    else:    
        return sigma * (2 - 5 * n)

def Dn(c):
    return alpha * c  + beta

# Primary Model
def Model(dx, dt = 1 / 12, time_iter_count = 13, X = 15): # dt is given in paper (under figure 1), time_iter_count = number of time steps
    
    x_intervals = int(X/dx) + 1

    fig, axes = plt.subplots(1, 2)

    # columns indicate spaces steps, rows indicate time steps
    n_data = np.zeros((13, x_intervals)) # First x_position is always 1; represents conentration outside of wound (position - 1)
    c_data = np.zeros((13, x_intervals)) # First x_position is always 1; represents conentration outside of wound (position - 1)
    
    # time loop
    for j in range(0, time_iter_count):
            
        prev_time = j - 1
        
        # Initialize initial values
        if j == 0:
            
            # MIGHT NEED TO CHANGE THIS BECAUSE IT'S SETTING 1 EQUAL TO THE LEFT BOUNDARY 
            n_data[0,:] = 1 
            c_data[0,:] = 1

        else:
            
            # starting the space loop
            for i in range(x_intervals):
                
                # running through the space steps except for the boundary
                if not (i == 0 or i == x_intervals-1):
                # Difference equations:
                    
                    # Cell Density
                    # MAYBE CHANGE TO i and i-1 instead of i+1 and i.
                    cell_migration = (alpha/dx**2) * (c_data[j-1, i+1] - c_data[j-1, i]) * (n_data[j-1, i+1] - n_data[j-1, i]) + (alpha * c_data[j-1, i] + beta) * (n_data[j-1, i+1] - 2 * n_data[j-1, i] + n_data[j-1, i-1]) / dx**2
                    mitotic_generation = (0.915 * c_data[j-1, i] + 0.0851) * n_data[j-1, i] * (2 - n_data[j-1, i])
                    natural_loss = n_data[j-1, i]
                    
                    n_data[j, i] = n_data[j-1, i] + dt * (cell_migration + mitotic_generation - natural_loss)
                    
                    # EGF Concentration
                    diffusion = Dc * (c_data[j-1, i+1] - 2*c_data[j-1, i] + c_data[j-1, i-1]) / dx**2
                    cell_production = A + B(n_data[j-1, i])
                    EGF_decay = (mu*n_data[j-1, i]*c_data[j-1, i]) / (c_hat + c_data[j-1, i]) + delta * c_data[j-1, i]
                    
                    c_data[j, i] = c_data[j-1, i] + dt * (diffusion + cell_production - EGF_decay)

                # # Experiment with removing these. It's possible that the model naturally approaches these values
                elif i == x_intervals - 1:
                    n_data[j, i] = 0
                    c_data[j, i] = 50
                    
    # Print data to file
    with open("output.txt", 'w') as f:
        
        for row in n_data:
            f.write(' '.join(map(str,row)) + '\n')
            
        f.write('\n')
        
        for row in c_data:
            f.write(' '.join(map(str,row)) + '\n')
            
    # Plot data
    for j in range(len(n_data)):
        axes[0].plot(np.linspace(0, 10, x_intervals), n_data[j], label = "Time = " + str(dt * j))
        axes[1].plot(np.linspace(0, 10, x_intervals), c_data[j], label = "Time = " + str(dt * j))

    axes[0].set_xlabel("X")
    axes[1].set_xlabel("X")

    axes[0].set_ylabel("Cell Density")
    axes[1].set_ylabel("EGF Concentration")

    axes[1].legend()
    plt.show()

    print("Resulting data is output to Output.txt file")
        
Model(dx, dt)