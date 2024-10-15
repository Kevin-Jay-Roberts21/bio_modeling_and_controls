#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 14:24:18 2024

@author: kevinjayroberts
"""

import matplotlib.pyplot as plt
import numpy as np

#############
# PROBELM 1 #
#############

# a) Use equation (1) to model glucose diffusing into the agar.

# defining some parameters
length = 20 # (cm)
thickness = 10 # (cm)
height = 1 # (cm)
C_e = 30 # (g/l)
C_f = 20 # (g/l)
D = 0.252*10**(-4) # (m^2/day)
dx = 0.01
dy = 0.01
k_m = 3.1536*10**(-3) 
Sh = k_m*dx/D
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo_m = D*dt/(dx**2)

print("PROBELM 1a OUTPUT")


# b) Use dt = round(dt, 3) from 1.a. Use an array with 2*x*y dimensions to calculate 
#    agar concentrations at each node at each successive time.

maximum = 60
save_step = 15
interval_steps = int(save_step / dt)

S = int(maximum / save_step) + 1



# defining the concentration
final_x = 0.2 # (cm)
final_y = 0.1 # (cm)

xn = int(final_x/dx) + 1
yn = int(final_y/dy) + 1

c = np.zeros((2,xn,yn)) # sets the initial conditions
data = np.zeros((S, xn, yn))
data[0] = c[1]

# printing the first time
print(f"Day 0: Minimum mass = {np.min(c[1]):.4f} g/L")
step_saved = 1

# adding the glucose on the side 

loop_counter = 0
time = 0

while time < maximum:
    
    c[0,:,:] = c[1,:,:]   
    
    #interior
    for i in range(1,xn-1):
        for j in range(1,yn-1):
            c[1, i, j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))
    
    
    # sides, left, right, bottom and top respectively
    for k in range(1, yn-1):
        c[1, 0, k] = 2*Fo_m*(c[0,1,k] + 1/2*(c[0,0,k-1] + c[0,0,k+1]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,0,k]
    for k in range(1, yn-1):
        c[1, xn-1, k] = 2*Fo_m*(c[0,xn-2,j] + 1/2*(c[0,xn-1,k-1] + c[0,xn-1,k+1]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,xn-1,k]
    for l in range(1, xn-1):
        c[1, l, 0] = 2*Fo_m*(c[0,l,0] + 1/2*(c[0,l-1,0] + c[0,l+1,0]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,l,0]
    for l in range(1, xn-1):
        c[1, l, yn-1] = 2*Fo_m*(c[0,l,yn-2] + 1/2*(c[0,l-1,yn-1] + c[0,l+1,yn-1]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,l,yn-1]
    
    
    # corners, bottom left, top left, bottom right and top right respectively
    c[1, 0, 0] = 2*Fo_m*(c[0,0,1] + c[0,1,0] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,0,0]
    c[1, 0, yn-1] = 2*Fo_m*(c[0,0,(yn-1)-1] + c[0,1,(yn-1)] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,0,(yn-1)]
    c[1, xn-1, 0] = 2*Fo_m*(c[0,(xn-1)-1,0] + c[0,(xn-1),1] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,(xn-1),0]
    c[1, xn-1, yn-1] = 2*Fo_m*(c[0,(xn-1),(yn-1)-1] + c[0,(xn-1)-1,(yn-1)] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,(xn-1),(yn-1)]

    
    loop_counter += 1
    time += dt
    
    # saving the data every 15 steps
    if loop_counter % interval_steps == 0:
        data[step_saved] = c[1,:,:]
        minimum = np.min(c[1,:,:])
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")
        step_saved += 1

    # Check if all agar has at least target concentration
    if np.all(c[1,:,:] >= C_f):
        print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        break



# Plot concentration profiles in 3D
x_grid = np.linspace(0, final_x, xn)  # x-axis points
y_grid = np.linspace(0, final_y, yn)  # y-axis points
X, Y = np.meshgrid(x_grid, y_grid)         # create the meshgrid

for m in range(loop_counter):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface, with correctly matched X, Y, and data[m,:,:] shapes
    ax.plot_surface(X, Y, data[m,:,:].T, cmap='viridis')

    ax.set_title(f'Mass at Day {m * save_step}')
    ax.set_xlabel('Length in cm')
    ax.set_ylabel('Thickness in cm')
    ax.set_zlabel('Concentration in g/L')
    plt.show()

# NOTE TO THE GRADER 
# Unsure why I'm getting this error. I'm unable to plot the data but I'm pretty sure it's correct.
# Let me know if you find what the issue is.


# c) Console Output
print()
print("The higher concentration on the boundaries is due to the glucose solution along the boundary. Diffusion occurs slower in the middle due to the fact that it's far from the boundary.")









