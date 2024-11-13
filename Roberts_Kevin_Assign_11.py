#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:27:15 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#############
# PROBLEM 1 #
#############

T_f = 15 # (C) # target temp
save_step = 10 # (save every 10 seconds)

# other constant params
p = 800 # (kg/m^3) desnity
C_p = 2400 # (J/(kg*K)) heat capacity
kh = 0.4 # (W/(m*k)) thermal conductivity
h = 100 # (W/(K*m^2)) heat transfer coefficient
T_e = 90 # (C) external temperature
alpha = kh/(p*C_p) # thermal diffusivity

dx = 0.001 # (m)
dy = 0.001 # (m)
dz = 0.001 # (m)

x_total = 0.010 # (m) 
y_total = 0.015 # (m)
z_total = 0.035 # (m)

Sh = round(k_m*dx/D,3) # sherwood number
Bi = h*dx/kh # biot number

interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))

dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)

dt = round(dt_max, 3)

Fo_m = D*dt/(dx**2)
Fourier_number = (alpha*dt)/dx**2
V = x_total*y_total*z_total
A = x_total*y_total

estimate_cooling_time = (p*C_p*V)*18/(h*A*(90-80))

t_x_1D = x_total**2/alpha
t_y_1D = y_total**2/alpha
t_z_1D = z_total**2/alpha
t_x_2D = x_total**2/(4*alpha)
t_y_2D = y_total**2/(4*alpha)
t_z_2D = z_total**2/(4*alpha)

# a) Console Output
print("PROBELM 1a OUTPUT")
print("x: 0.010m")
print("y: 0.015m")
print("Aspect Ratio: 10:15 = 2:3")
print("Thermal Diffusivity: " + str(alpha))
print("Biot number: " + str(Bi))
print("dt_max: " + str(dt_max))
print("Fourier number: " + str(Fourier_number))
print("Final Cooling Time: " + str())
print("Estimates of the heat transfer time in x, y and z dimensions:")
print("t_x_1D = " + str(t_x_1D))
print("t_y_1D = " + str(t_y_1D))
print("t_z_1D = " + str(t_z_1D))
print("Estimate of the heat transfer in two dimensions for x, y and z:")
print("t_x_2D = " + str(t_x_1D))
print("t_y_2D = " + str(t_y_2D))
print("t_z_2D = " + str(t_z_2D))

# b) Figure 1.a
def run_2D_heat_transfer(x, y, z):
    pass




























def run_3D_heat_transfer(x, y, z):
    pass

def plot_min_heats(time_grid, dimless_time_grid, mins, cube_id):
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    # First subplot for the original time plot
    ax1.plot(time_grid, mins, marker='o', label='Cube ' + cube_id + ' Heat minimums')

    # Customizing the first subplot
    ax1.set_title('Heat Cube ' + cube_id + ' Mins at every 10 seconds')
    ax1.set_xlabel("Time (seconds)")
    ax1.set_ylabel('Minimum Heat in Cube ' + cube_id + ' Values')
    ax1.set_xticks(time_grid)  # Set x-ticks to match the time values
    ax1.legend()
    ax1.grid()

    # Second subplot for the dimensionless time plot
    ax2.plot(dimless_time_grid, mins, marker='o', label='Cube ' + cube_id + ' Heat minimums')

    # Customizing the second subplot
    ax2.set_title('Heat Cube ' + cube_id + ' Mins (Dimensionless Time)')
    ax2.set_xlabel("Dimensionless Time")
    ax2.set_ylabel('Minimum Heat in Cube ' + cube_id + ' Values')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_xticks(dimless_time_grid)  # Set x-ticks to match the time values
    ax2.legend()
    ax2.grid()

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

def display_heat_data(cube_data, cube_lengths, cube_id):
    print("#################")
    print("Results of Cube " + cube_id)
    print("x: " + cube_lengths + " cm")
    print("y: " + cube_lengths + " cm")
    print("z: " + cube_lengths + " cm")
    print("Thermal diffusivity: " + str(cube_data[4])) 
    print("Biot Number: " + str(round(cube_data[3], 3)))
    
    # printing the minimum concentration at a certain position for every 10 seconds
    for i in range(cube_data[6]):
        
        if i == cube_data[6]-1:
            time_of_min_heat = str(cube_data[2])
        else:
            time_of_min_heat = str(i*10)
        
        position = np.unravel_index(np.argmin(cube_data[0][i]), cube_data[0][i].shape)
        print("T = " + time_of_min_heat + " seconds:   min heat: " + str(round(np.min(cube_data[0][i]), 3)) + "    position: " + str(position))

cube1_data = run_heat_transfer(cube1_width, cube1_height, cube1_thickness)
cube2_data = run_heat_transfer(cube2_width, cube2_height, cube2_thickness)
cube3_data = run_heat_transfer(cube3_width, cube3_height, cube3_thickness)
