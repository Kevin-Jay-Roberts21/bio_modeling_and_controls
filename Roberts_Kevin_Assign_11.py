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
T_i = 33 # (C) # initial temp 
save_step = 10 # (save every 10 seconds)

# other constant params
p = 800 # (kg/m^3) desnity
C_p = 2400 # (J/(kg*K)) heat capacity
kh = 0.4 # (W/(m*k)) thermal conductivity
h = 100 # (W/(K*m^2)) heat transfer coefficient
T_e = 90 # (C) external temperature
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
print("Final Cooling Time: " + str(round(estimate_cooling_time)))
print("Estimates of the heat transfer time in x, y and z dimensions:")
print("t_x_1D = " + str(round(t_x_1D)))
print("t_y_1D = " + str(round(t_y_1D)))
print("t_z_1D = " + str(round(t_z_1D)))
print("Estimate of the heat transfer in two dimensions for x, y and z:")
print("t_x_2D = " + str(round(t_x_1D)))
print("t_y_2D = " + str(round(t_y_2D)))
print("t_z_2D = " + str(round(t_z_2D,3)))
print()
print()

# b) Figure 1.a

xn = int(x_total/dx)
yn = int(y_total/dy)

c = np.zeros((2, xn+1, yn+1)) # sets the initial conditions
data = np.zeros((1, xn+1, yn+1))
data[0] = c[0]

step_saved = 0

loop_counter = 0
time = 0
end_time = 0
minimum = 0

while minimum < T_f:
    
    c[0,:,:] = c[1,:,:]   
    
    for i in range(0,xn+1):
        for j in range(0,yn+1):
            if i == 0 and j == 0:  # upper left corner
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 0.5*2*Sh*T_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0 and j == yn:  # upper right corner
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 0.5*2*Sh*T_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == 0:  # bottom left corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 0.5*2*Sh*T_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == yn:  # bottom right corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 0.5*2*Sh*T_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0:  # top
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + 0*Sh*T_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif j == 0:  # left side
                c[1,i,j] = 2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + 0*Sh*T_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + 0*Sh*T_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
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
        minimum = np.min(c[1,:,:])
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")

    # Check if all agar has at least target concentration
    if np.all(c[1,:,:] >= T_f):
        # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        end_time = dt*loop_counter
        break

for i in range(len(data)-1):
    minimum = np.min(data[i,:,:])
    time = i*save_step
    print(f"{int(time)} secs: Snickers bar temp = {minimum:.3f} C")
end_min = np.min(data[len(data)-1,:,:])
print(f"{int(end_time)} secs: Snickers bar temp = {end_min:.3f} C")


# c) Explain
print("PROBELM 1c OUTPUT")
print("i) ")
print("ii) ")
print()
print()


# d) Figure 1b

# e) Explain
print("PROBLEM 1e OUTPUT")
print("i) ")
print("ii) ")
print("iii) ")
print()
print()

#############
# PROBLEM 2 #
#############

# a) Console output
print("PROBLEM 2a OUTPUT")
print("alpha: " + str())
print("Biot number: " + str())
print("T_max: " + str())
# print coordinates for ...
# print perc diff
print()
print()

# b) Explain
print("PROBLEM 2b OUTPUT")
print("i) ")
print("ii) ")
print("iii) ")

def run_3D_heat_transfer(x, y, z):
    pass


