#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 14:28:53 2024

@author: kevinjayroberts
"""

import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(precision=3, edgeitems=3, suppress=True)

#############
# PROBLEM 1 #
#############

# a) console output
def get_agar_data(total_x):

    # defining variables
    C_e = 30 # (g/l)
    C_f = 20 # (g/l)
    D = 0.252*(10**(-4)) # (m^2/day)
    k_m = 3.153*10**(-3) # (m/day)
    dx = 0.01 # (m)
    dy = 0.01 # (m)
    # total_x is a function parameter
    total_y = 0.1 # (m)
    
    Sh = k_m*dx/D
    interior_dt = dx**2/(4*D)
    side_dt = dx**2/(2*D*(2 + Sh))
    corner_dt = dx**2/(4*D*(1 + Sh))
    dt_maxes = [interior_dt, side_dt, corner_dt] 
    dt_max = min(dt_maxes)
    dt = round(dt_max, 3)
    Fo_m = D*dt/(dx**2)
    
    maximum = 90
    save_step = 15
    interval_steps = int(save_step / dt)
    S = int(maximum / save_step) + 1
    
    xn = int(total_x/dx) + 1
    yn = int(total_y/dy) + 1
    
    c = np.zeros((2,xn,yn)) # sets the initial conditions
    data = np.zeros((S, xn, yn))
    data[0] = c[1]
    
    # printing the first time
    #print(f"Day 0: Minimum mass = {np.min(c[1]):.4f} g/L")
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
            data[step_saved] = c[1]
            minimum = np.min(c[1])
            #print(f"Day {int(time)}: Minimum mass concentration = {minimum:.4f} g/L")
            step_saved += 1
        
        i# Check if all agar has at least target concentration
        if np.all(c[1,:,:] >= C_f):
            #print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
            break
        
    
    information = [c, data, time, loop_counter]
    
    return information

# data for different p's
diff_p_data = []

p = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30] # all in (m)

print("PROBELM 1a OUTPUT:")
for i in range(len(p)):
    
    x_length = 0.2+p[i]
    p_data = get_agar_data(x_length)
    
    diff_p_data.append(p_data)
    print("x[m]: " + str(round(p_data[0][:,1][1][len(p_data[0][:,1][1])-1],3)) + ", aspect ratio: " + str(int(x_length*100)) + "x10x1," + " final time: " + str(round(p_data[2],3)) + ", perc diff: " + str(1))





# b) Figure 1

# c) Explain
print("i) ")
print("ii) ")
print("iii) ")
print("iv) ")
print("v) ")
print("vi) ")
print("vii) ")
print("viii) ")


#############
# PROBLEM 2 #
#############

# a) console output

# b) Figure 2.a

# c) Explain
print("i) ")
print("ii) ")
print("iii) ")
print("iv) ")

# d) Figure 2.b

# e) Explain
print("i) ")
print("ii) ")
print("iii) ")


#############
# PROBLEM 3 #
#############

