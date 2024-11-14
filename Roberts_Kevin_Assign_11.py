#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:27:15 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def perc_diff(A, B):
    return (np.abs(A - B)*100)/((A+B)/2)

#############
# PROBLEM 1 #
#############

T_f = 15 # (C) # target temp
T_i = 33 # (C) # initial temp 

save_step = 10 # (save every 10 seconds)

# other constant params
p = 800 # (kg/m^3) desnity
C_p = 2400 # (J/(kg*K)) heat capacity
kh = 0.4 # (W/(m*k)) thermal conductivity
h = 100 # (W/(K*m^2)) heat transfer coefficient
T_e = 5 # (C) external temperature
alpha = kh/(p*C_p) # thermal diffusivity

D = alpha

dx = 0.001 # (m)
dy = 0.001 # (m)
dz = 0.001 # (m)

x_total = 0.010 # (m) 
y_total = 0.015 # (m)
z_total = 0.035 # (m)

Bi = h*dx/kh # biot number
Sh = Bi # sherwood number

interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))

dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)

dt = round(dt_max,3)

Fo = D*dt/(dx**2)
V = x_total*y_total*z_total
A = x_total*y_total

estimate_cooling_time = (p*C_p*V)*18/(h*A*(90-80))

L_x = 10 * 0.001
L_y = 15/2 * 0.001
L_z = 35/2 * 0.001

t_x_1D = L_x**2/alpha
t_y_1D = L_y**2/alpha
t_z_1D = L_z**2/alpha
t_x_2D = L_x**2/(4*alpha)
t_y_2D = L_y**2/(4*alpha)
t_z_2D = L_z**2/(4*alpha)

save_step = 50
interval_steps = int(save_step / dt)

# a) Console Output
print("PROBELM 1a OUTPUT")
print("x: 0.010m")
print("y: 0.015m")
print("Aspect Ratio: 10:15 = 2:3")
print("Thermal Diffusivity: " + str(alpha))
print("Biot number: " + str(Bi))
print("dt_max: " + str(dt_max))
print("Fourier number: " + str(Fo))
print("Final Cooling Time: " + str(round(estimate_cooling_time,3)))
print("Estimates of the heat transfer time in x, y and z dimensions:")
print("t_x_1D = " + str(round(t_x_1D,3)))
print("t_y_1D = " + str(round(t_y_1D,3)))
print("t_z_1D = " + str(round(t_z_1D,3)))
print("Estimate of the heat transfer in two dimensions for x, y and z:")
print("t_x_2D = " + str(round(t_x_2D,3)))
print("t_y_2D = " + str(round(t_y_2D,3)))
print("t_z_2D = " + str(round(t_z_2D,3)))
print()
print()

#################
# b) Figure 1.a #
#################
xn = int(x_total/dx)
yn = int(y_total/dy)

c = T_i * np.ones((2, xn+1, yn+1)) # sets the initial conditions
data = np.zeros((1, xn+1, yn+1))
data[0] = c[0]

step_saved = 0

loop_counter = 0
time = 0
end_time = 0

while True:
    
    c[0,:,:] = c[1,:,:]   
    
    for i in range(0,xn+1):
        for j in range(0,yn+1):
            if i == 0 and j == 0:  # upper left corner
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 2*Sh*T_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,i,j]
            elif i == 0 and j == yn:  # upper right corner
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 2*Sh*T_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == 0:  # bottom left corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 0.5*2*Sh*T_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == yn:  # bottom right corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 0.5*2*Sh*T_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0:  # top
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*T_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == 0:  # left side
                c[1,i,j] = 2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*T_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*T_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif i == xn:  # bottom
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + 0*Sh*T_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            else:
                c[1,i,j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))

    loop_counter += 1
    time += dt
    
    # saving the data every 50 steps
    if loop_counter % interval_steps == 0:
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")

    # Check if all temp cools to 15 C or lower
    if np.all(c[1,:,:] <= T_f):
        # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        end_time = dt*loop_counter
        break

# collecting all of the mins and maxs
mins = []
maxs = []

# Compute minimum and maximum values for each time step
for i in range(len(data)):
    mins.append(np.min(data[i]))
    maxs.append(np.max(data[i]))
print("Temp mins:")
# Print minimum values at each time step
for i in range(len(mins) - 1):
    minimum = mins[i]
    time = i * save_step
    print(f"{int(time)} secs: Snickers bar temp (min) = {minimum:.3f} C")
end_min = np.min(data[len(data) - 1, :, :])
print(f"{int(end_time)} secs: Snickers bar temp (min) = {end_min:.3f} C")
print()

print("Temp maxes:")
# Print maximum values at each time step
for i in range(len(maxs) - 1):
    maximum = maxs[i]
    time = i * save_step
    print(f"{int(time)} secs: Snickers bar temp (max) = {maximum:.3f} C")
end_max = np.max(data[len(data) - 1, :, :])
print(f"{int(end_time)} secs: Snickers bar temp (max) = {end_max:.3f} C")
print()
print()

# Plotting parameters
L = L_y  # (m)
t = np.arange(0, (len(data) - 1) * save_step, save_step)
t = np.append(t, end_time)
dim_less_time = [round(i * (D / L**2), 3) for i in t]

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# First subplot for the original time plot
ax1.plot(t, mins, marker='o', label='Temperature minimums')
ax1.plot(t, maxs, marker='o', label='Temperature maximums')  # Add the max plot

# Customizing the first subplot
ax1.set_title("Plotting Temp Mins and Maxs at every 50 seconds")
ax1.set_xlabel("Time (seconds)")
ax1.set_ylabel("Temp values (C)")
ax1.set_xticks(t)  # Set x-ticks to match the time values
ax1.legend()
ax1.grid()

# Second subplot for the dimensionless time plot
ax2.plot(dim_less_time, mins, marker='o', label='Temperature minimums')
ax2.plot(dim_less_time, maxs, marker='o', label='Temperature maximums')  # Add the max plot

# Customizing the second subplot
ax2.set_title("Plotting Temp Mins and Maxs (Dimensionless Time)")
ax2.set_xlabel("Dimensionless Time")
ax2.set_ylabel("Temp values (C)")
ax2.tick_params(axis='x', rotation=45)
ax2.set_xticks(dim_less_time)  # Set x-ticks to match the time values
ax2.legend()
ax2.grid()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()


##############
# c) Explain #
##############
print("PROBELM 1c OUTPUT")
print("i) The difference in max temp of the snickers bar differs very much from the minimum initially, however as time increases, they tend to get closer to each other.")
print("ii) The difference in the estimated temp and the actual temp for when the snickers bar reaches 15 C is not very big. This is as expected however, we now wish to compare with the 3D model.")
print()
print()


# d) Figure 1b
x_grid = np.linspace(0, x_total, xn+1)  # x-axis points
y_grid = np.linspace(0, y_total, yn+1)  # y-axis points
X, Y = np.meshgrid(x_grid, y_grid)      # create the meshgrid

for m in range(len(data)-1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface, with correctly matched X, Y, and data[m,:,:] shapes
    ax.plot_surface(X, Y, data[m,:,:].T, cmap='viridis')

    ax.set_title(f'Temp at Day {m * save_step}')
    ax.set_xlabel('Length in cm')
    ax.set_ylabel('Thickness in cm')
    ax.set_zlabel('Temp in C')
    plt.show()
# plotting the last day
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface, with correctly matched X, Y, and data[m,:,:] shapes
ax.plot_surface(X, Y, data[len(data)-1,:,:].T, cmap='viridis')

ax.set_title(f'Temp at Day {int(end_time)}')
ax.set_xlabel('Length in cm')
ax.set_ylabel('Thickness in cm')
ax.set_zlabel('Temperature in C')
plt.show()


# e) Explain
print("PROBLEM 1e OUTPUT")
print("i) The rapid initial temperature change is followed by slower convergence. The changes to the boundary conditions greatly affect the heat transfer rates.")
print("ii) The one-dimensional model may be sufficient for thin bars with uniform corss-sections and fast heat transfer times.")
print("iii) Increase the convective heat transfer coefficient, thermal conductivity, or the surface area-to-volume ratio.")
print()
print()

#############
# PROBLEM 2 #
#############

# resetting all of the values to ensure correctness

T_f = 15 # (C) # target temp
T_i = 33 # (C) # initial temp 

save_step = 10 # (save every 10 seconds)

# other constant params
p = 800 # (kg/m^3) desnity
C_p = 2400 # (J/(kg*K)) heat capacity
kh = 0.4 # (W/(m*k)) thermal conductivity
h = 100 # (W/(K*m^2)) heat transfer coefficient
Te = 5 # (C) external temperature
alpha = kh/(p*C_p) # thermal diffusivity

D = alpha

dx = 0.001 # (m)
dy = 0.001 # (m)
dz = 0.001 # (m)

x_total = 0.010 # (m) 
y_total = 0.015 # (m)
z_total = 0.035 # (m)

Bi = h*dx/kh # biot number
Sh = Bi # sherwood number

interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))

dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)

dt = round(dt_max,3)
dt = 0.5
Fo = D*dt/(dx**2)
V = x_total*y_total*z_total
A = x_total*y_total

estimate_cooling_time = (p*C_p*V)*18/(h*A*(90-80))

L_x = 10 * 0.001
L_y = 15/2 * 0.001
L_z = 35/2 * 0.001

t_x_1D = L_x**2/alpha
t_y_1D = L_y**2/alpha
t_z_1D = L_z**2/alpha
t_x_2D = L_x**2/(4*alpha)
t_y_2D = L_y**2/(4*alpha)
t_z_2D = L_z**2/(4*alpha)

save_step = 50
interval_steps = int(save_step / dt)

xn = int(x_total/dx)
yn = int(y_total/dy)
zn = int(z_total/dz)

c = T_i * np.ones((2, zn+1, yn+1, xn+1))  # sets the initial conditions
data = np.zeros((1, zn+1, yn+1, xn+1))
data[0] = c[0]

step_saved = 0

loop_counter = 0
time = 0
end_time = 0
while True: 
    c[0,:,:,:] = c[1,:,:,:]
    for i in range(0, xn+1):
        for j in range(0, yn+1):
            for k in range(0, zn+1):
                ###########
                # CORNERS #
                ###########
                if i == 0 and j == 0 and k == 0: # left, top, forward
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == 0 and j == 0 and k == zn: # left, top, back
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == xn and j == 0 and k == zn: # right, top, back
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif i == xn and j == 0 and k == 0: # right, top, forward
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == 0 and j == yn and k == 0: # left, bottom, forward
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == 0 and j == yn and k == zn: # left, bottom, back
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == xn and j == yn and k == zn: # right, bottom, back
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif i == xn and j == yn and k == 0: # right, bottom, forward
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                #########
                # EDGES #                   
                #########
                elif i == 0 and k == 0: # left, forward edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == 0 and k == zn: # left, back edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif i == xn and k == zn: # right, back edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif i == xn and k == 0: # right, forward edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif j == 0 and k == 0: # front, upper edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif j == 0 and k == zn: # back, upper edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif j == yn and k == zn: # back, lower edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif j == yn and k == 0: # front, lower edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))
                elif i == 0 and j == 0: # left, upper edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == xn and j == 0: # right, upper edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif i == xn and j == yn: # right, lower edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif i == 0 and j == yn: # left, lower edge
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                #########
                # FACES #
                #########                  
                elif k == 0: # front face
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0)/ dz**2 +
                        (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                        (c[0,k,j,i+1] - 2 * c[0,k,j,i] +c[0,k,j,i-1])/ dx**2) +
                            2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif k == zn: # back face
                     c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                         (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0)/ dz**2 +
                         (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                         (c[0,k,j,i+1] - 2 * c[0,k,j,i] +c[0,k,j,i-1])/ dx**2) +
                             2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                             0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                             0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                                      
                elif i == 0: # left face
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                        (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                        (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0)/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif i == xn: # right face
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                        (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                        (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0)/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                elif j == 0: # top face
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                        (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0)/ dy**2 +
                        (c[0,k,j,i+1] - 2 * c[0,k,j,i] +c[0,k,j,i-1])/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                elif j == yn: # bottom face
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                        (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0)/ dy**2 +
                        (c[0,k,j,i+1] - 2 * c[0,k,j,i] +c[0,k,j,i-1])/ dx**2) +
                            0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                            0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                # CENTER                            
                else:
                    c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                        (c[0,k+1,j,i] - 2 * c[0,k,j,i] + c[0,k-1,j,i])/ dz**2 +
                        (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                        (c[0,k,j,i+1] - 2 * c[0,k,j,i] + c[0,k,j,i-1])/ dx**2) )
   
    loop_counter += 1
    
    if loop_counter % interval_steps == 0:
        step_saved += 1
        data = np.append(data, np.zeros((1, zn+1, yn+1, xn+1)), axis=0)
        data[step_saved] = c[1,:,:,:]
       
    # print(np.max(c[1]))
       
    if np.all(c[1,:,:,:] <= T_f):
        step_saved += 1
        data = np.append(data, np.zeros((1, zn+1, yn+1, xn+1)), axis=0)
        data[step_saved] = c[1,:,:,:]
        end_time1 = loop_counter * dt
        break

    if loop_counter > 500: # assuming heat did not decay
        break

end_time1 = dt*loop_counter

perc_difference = perc_diff(end_time, end_time1)

# a) Console output
print("PROBLEM 2a OUTPUT")
print("alpha: " + str(alpha))
print("Biot number: " + str(Bi))

mins = []
maxs = []

# Compute minimum and maximum values for each time step
for i in range(len(data)-1):
    mins.append(np.min(data[i]))
    maxs.append(np.max(data[i]))

# printing the max temp at a certain position for every 50 seconds
for i in range(len(data)-1):
    
    if i == len(data)-1:
        time_of_max_heat = str(end_time)
    else:
        time_of_max_heat = str(i*50)
    
    position = np.unravel_index(np.argmax(data[i]), data[i].shape)
    print("T = " + time_of_max_heat + " seconds:   max heat: " + str(round(np.max(data[i]), 3)) + "    position: " + str(position))

print("Percentage Difference: " + str(perc_difference))
print()
print()


# Plotting parameters
L = L_y  # (m)
t = np.arange(0, (len(data) - 1) * save_step, save_step)
dim_less_time = [round(i * (D / L**2), 3) for i in t]

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# First subplot for the original time plot
ax1.plot(t, mins, marker='o', label='Temperature minimums')
ax1.plot(t, maxs, marker='o', label='Temperature maximums')  # Add the max plot

# Customizing the first subplot
ax1.set_title("Plotting 3D Temp Mins and Maxs at every 50 seconds")
ax1.set_xlabel("Time (seconds)")
ax1.set_ylabel("Temp values (C)")
ax1.set_xticks(t)  # Set x-ticks to match the time values
ax1.legend()
ax1.grid()

# Second subplot for the dimensionless time plot
ax2.plot(dim_less_time, mins, marker='o', label='Temperature minimums')
ax2.plot(dim_less_time, maxs, marker='o', label='Temperature maximums')  # Add the max plot

# Customizing the second subplot
ax2.set_title("Plotting 3D Temp Mins and Maxs (Dimensionless Time)")
ax2.set_xlabel("Dimensionless Time")
ax2.set_ylabel("Temp values (C)")
ax2.tick_params(axis='x', rotation=45)
ax2.set_xticks(dim_less_time)  # Set x-ticks to match the time values
ax2.legend()
ax2.grid()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()

# b) Explain
print("PROBLEM 2b OUTPUT")
print("i) The heat in the 3D model cools faster than the heat in the 2D model.")
print("ii) The 3D model takes much langer to run due to the many if statements in the loop. Additionally, for the 2D model, we don't even have a 3rd loop so it runs much faster.")
print("iii) It may be more useful to use a 3D model for a short time, however for long times, it may be better to use the 2D model to avoid computational expense.")



