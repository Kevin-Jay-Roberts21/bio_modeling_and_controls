#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:28:48 2024

@author: kevinjayroberts
"""
import numpy as np
import matplotlib.pyplot as plt

#############
# PROBLEM 1 #
#############

# a) Console Output

z = 0.4 # length (m)
y = 0.1 # thickness (m)
x = 0.1 # height (m)

D = 0.252*(10**-4) # (m^2/day)
dx = 0.01 # (cm)
dy = 0.01 # (cm)
k_m = 3.1536*10**(-3) 
C_e = 20 # (g/l)
C_f = 15 # (g/l)
Sh = k_m*dx/D
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo_m = D*dt/(dx**2)

maximum = 60
save_step = 15
interval_steps = int(save_step / dt)

S = int(maximum / save_step) + 1

# defining the concentration
final_x = 0.1 # (m)
final_y = 0.1 # (m)

xn = int(final_x/dx) + 1
yn = int(final_y/dy) + 1

c = np.zeros((2,xn,yn)) # sets the initial conditions
data = np.zeros((S, xn, yn))
data[0] = c[1]

# printing the first time
print(f"Day 0: Minimum analgesic = {np.min(c[1]):.4f} g/L")
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
        c[1, xn-1, k] = 2*Fo_m*(c[0,xn-2,j] + 1/2*(c[0,xn-1,k-1] + c[0,xn-1,k+1]) + 0*Sh*C_e) + (1 - 4*Fo_m - 0*2*Fo_m*Sh)*c[0,xn-1,k]
    for l in range(1, xn-1):
        c[1, l, 0] = 2*Fo_m*(c[0,l,0] + 1/2*(c[0,l-1,0] + c[0,l+1,0]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,l,0]
    for l in range(1, xn-1):
        c[1, l, yn-1] = 2*Fo_m*(c[0,l,yn-2] + 1/2*(c[0,l-1,yn-1] + c[0,l+1,yn-1]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,l,yn-1]
    
    # corners, bottom left, top left, bottom right and top right respectively
    c[1, 0, 0] = 2*Fo_m*(c[0,0,1] + c[0,1,0] + 0.5*2*Sh*C_e) + (1 - 4*Fo_m - 0.5*4*Fo_m*Sh)*c[0,0,0]
    c[1, 0, yn-1] = 2*Fo_m*(c[0,0,(yn-1)-1] + c[0,1,(yn-1)] + 0.5*2*Sh*C_e) + (1 - 4*Fo_m - 0.5*4*Fo_m*Sh)*c[0,0,(yn-1)]
    c[1, xn-1, 0] = 2*Fo_m*(c[0,(xn-1)-1,0] + c[0,(xn-1),1] + 0.5*2*Sh*C_e) + (1 - 4*Fo_m - 0.5*4*Fo_m*Sh)*c[0,(xn-1),0]
    c[1, xn-1, yn-1] = 2*Fo_m*(c[0,(xn-1),(yn-1)-1] + c[0,(xn-1)-1,(yn-1)] + 0.5*2*Sh*C_e) + (1 - 4*Fo_m - 0.5*4*Fo_m*Sh)*c[0,(xn-1),(yn-1)]

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
        # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        break


# Plot concentration profiles in 3D
x_grid = np.linspace(0, final_x, xn)  # x-axis points
y_grid = np.linspace(0, final_y, yn)  # y-axis points
X, Y = np.meshgrid(x_grid, y_grid)         # create the meshgrid

for m in range(S):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface, with correctly matched X, Y, and data[m,:,:] shapes
    ax.plot_surface(X, Y, data[m,:,:].T, cmap='viridis')

    ax.set_title(f'Mass at Day {m * save_step}')
    ax.set_xlabel('Length in cm')
    ax.set_ylabel('Thickness in cm')
    ax.set_zlabel('Concentration in g/L')
    plt.show()


# printing all of the output stuff
print()
print("The higher concentration on the boundaries is due to the glucose solution along the boundary. Diffusion occurs slower in the middle due to the fact that it's far from the boundary.")

print("PROBELM 1a OUTPUT:")
for i in range(len(p)):
    
    x_length = 0.2+p[i]
    p_data = get_agar_data(x_length)
    
    diff_p_data.append(p_data)
    
    perc_diff_value = round(perc_diff(prev_data[2], p_data[2]), 3)
    transfer_time = round(p_data[2],3)
    
    
    if i == 0:
        print("x[m]: " + str(round(c[0][:,1][1][len(c[0][:,1][1])-1],3)) + ", aspect ratio: 40x10x10," + " final time: " + str(transfer_time))
    else:
        print("x[m]: " + str(round(c[0][:,1][1][len(c[0][:,1][1])-1],3)) + ", aspect ratio: 40x10x10," + " final time: " + str(transfer_time))
    prev_data = p_data
    perc_differences.append(perc_diff_value)
    transfer_times.append(transfer_time)



