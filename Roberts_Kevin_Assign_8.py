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
def perc_diff(A, B):
    return (np.abs(A - B)*100)/((A+B)/2)

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

p0 = get_agar_data(0.2)
prev_data = p0

perc_differences = []
transfer_times = []

print("PROBELM 1a OUTPUT:")
for i in range(len(p)):
    
    x_length = 0.2+p[i]
    p_data = get_agar_data(x_length)
    
    diff_p_data.append(p_data)
    
    perc_diff_value = round(perc_diff(prev_data[2], p_data[2]), 3)
    transfer_time = round(p_data[2],3)
    
    
    if i == 0:
        print("x[m]: " + str(round(p_data[0][:,1][1][len(p_data[0][:,1][1])-1],3)) + ", aspect ratio: " + str(int(x_length*100)) + "x10x1," + " final time: " + str(transfer_time) + ", perc diff: nan")
    else:
        print("x[m]: " + str(round(p_data[0][:,1][1][len(p_data[0][:,1][1])-1],3)) + ", aspect ratio: " + str(int(x_length*100)) + "x10x1," + " final time: " + str(transfer_time) + ", perc diff: " + str(perc_diff_value))
    prev_data = p_data
    perc_differences.append(perc_diff_value)
    transfer_times.append(transfer_time)


# b) Figure 1

x = np.array(p)

fig, ax1 = plt.subplots()
fig.suptitle("Max of x and y for each dt.")
line1 = ax1.scatter(x, perc_differences, c="blue", marker="D", label="perc difference") # Plot (x^3 vs x)
ax2 = ax1.twinx()
line2 = ax2.scatter(x, transfer_times, c="red", linestyle="solid", label="transfer time")
ax1.set_xlabel('Aspect ratio')
ax1.set_ylabel('percent difference', color = 'blue')
ax2.set_ylabel('transfer time', color = 'red')
ax1.tick_params('y', colors = 'blue')
ax2.tick_params('y', colors = 'red')
ax1.legend((line1,line2),('perc difference', 'transfer time'), loc = 'upper right', bbox_to_anchor=(1,-0.1))
plt.tight_layout()
plt.show()



# c) Explain
print()
print("PROBLEM 1c OUTPUT")
print("i) The mass transfer increase time increaseed but plateaued, and the percentage difference approched 0.")
print("ii) Yes, in very narrow environments. Mass transfer tends to have less of a diffusive behaviro and more of an advective behavior.")
print("iii) From the figure, we see the characteristic length is anything ax10x1 where a is anything 40 or higher.")
print("iv) A little under 63 days.")
print("v) As the aspect ratio increases, the characteristic time approach the same final time.")
print("vi) They should all be 1 or close to 1.")
print("vii) For infinite dimensions we expect so see no effect from the boundary condtions, and the opposite is true for controlled dimensions, where boundary conditions take effect.")
print("viii) Similar to what was said earilier: the more narrow the environment, the less dffusive behvior and vice versa.")
print()

#############
# PROBLEM 2 #
#############
print("PROBELM 2a OUTPUT")


# a) console output
tissue_slab_1 = [.4, .2, .01]
tissue_slab_2 = [.48, .12, .01]
tissue_slab_3 = [.56, .14, .01]

def get_mass_transfer_data(total_x, total_y):

    # defining variables
    C_e = 30 # (g/l)
    C_f = 20 # (g/l)
    D = 0.252*(10**(-4)) # (m^2/day)
    k_m = 3.153*10**(-3) # (m/day)
    dx = 0.01 # (m)
    dy = 0.01 # (m)
    
    Sh = k_m*dx/D
    interior_dt = dx**2/(4*D)
    side_dt = dx**2/(2*D*(2 + Sh))
    corner_dt = dx**2/(4*D*(1 + Sh))
    dt_maxes = [interior_dt, side_dt, corner_dt] 
    dt_max = min(dt_maxes)
    dt = round(dt_max, 3)
    Fo_m = D*dt/(dx**2)
    
    maximum = 230
    save_step = 15
    interval_steps = int(save_step / dt)
    S = int(maximum / save_step) + 1
    
    xn = int(total_x/dx) + 1
    yn = int(total_y/dy) + 1
    
    c = np.zeros((2,xn,yn)) # sets the initial conditions
    data = np.zeros((S, xn, yn))
    data[0] = c[1]
    
    # printing the first time
    print(f"Day 0: Minimum agar = {np.min(c[1]):.4f} g/L")
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
            print(f"Day {int(time)}: Minimum agar = {minimum:.4f} g/L")
            step_saved += 1
        
        # Check if all agar has at least target concentration
        if np.all(c[1,:,:] >= C_f):
            print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
            break
    
    information = [data, c[1][:,1][xn-1], c[1][:,1][0], str(int(total_x*100))+"x"+str(int(total_y*100))+"x1", time, Sh]
    
    return information

slab_1_data = get_mass_transfer_data(tissue_slab_1[0], tissue_slab_1[1])
print("x[m]: " + str(round(slab_1_data[1], 3)) + ", y[m]: " + str(round(slab_1_data[2], 3)) + ", aspect ratio: " + slab_1_data[3] + ", final time: " + str(round(slab_1_data[4],3)) + ", Sh: " + str(round(slab_1_data[5],3)))
print()

slab_2_data = get_mass_transfer_data(tissue_slab_2[0], tissue_slab_2[1])
print("x[m]: " + str(round(slab_2_data[1], 3)) + ", y[m]: " + str(round(slab_2_data[2], 3)) + ", aspect ratio: " + slab_2_data[3] + ", final time: " + str(round(slab_2_data[4],3)) + ", Sh: " + str(round(slab_2_data[5],3)))
print()

slab_3_data = get_mass_transfer_data(tissue_slab_3[0], tissue_slab_3[1])
print("x[m]: " + str(round(slab_3_data[1], 3)) + ", y[m]: " + str(round(slab_3_data[2], 3)) + ", aspect ratio: " + slab_3_data[3] + ", final time: " + str(round(slab_3_data[4],3)) + ", Sh: " + str(round(slab_3_data[5],3)))
print()


# b) Figure 2.a
# getting the minimums of each slab at each 15 days
slab_1_mins = []
slab_2_mins = []
slab_3_mins = []

# defining the length of each array
slab_1_data_len = len(slab_1_data[0])
slab_2_data_len = len(slab_2_data[0])
slab_3_data_len = len(slab_3_data[0])

max_time_of_all_slabs = max([slab_1_data_len, slab_2_data_len, slab_3_data_len])

for i in range(slab_1_data_len-1):
    
    if np.min(slab_1_data[0][i]) == 0 and i > 0:
        break
    else:
        slab_1_mins.append(np.min(slab_1_data[0][i]))
for i in range(slab_2_data_len-1):
    if np.min(slab_2_data[0][i]) == 0 and i > 0:
        break
    else:
        slab_2_mins.append(np.min(slab_2_data[0][i]))
for i in range(slab_3_data_len-1):
    if np.min(slab_3_data[0][i]) == 0 and i > 0:
        break
    else:
        slab_3_mins.append(np.min(slab_3_data[0][i]))


t = [0, 14, 29, 44, 59, 74, 89, 104, 119, 134, 149, 164, 179, 194, 209]

# Assuming slab_1_mins, slab_2_mins, slab_3_mins, and t are already defined

# Converting to dimensionless time
D = 0.252 * (10**(-4))  # (m^2/day)
L1 = 0.40
L2 = 0.48
L3 = 0.56

# Characteristic length
L = (L1 + L2 + L3) / 3

dim_less_time = [round(i * (D / L**2), 3) for i in t]

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# First subplot for the original time plot
ax1.plot(t, slab_1_mins, marker='o', label='Slab 1 Minimum')
ax1.plot(t[:len(slab_2_mins)], slab_2_mins, marker='s', label='Slab 2 Minimum')
ax1.plot(t[:len(slab_3_mins)], slab_3_mins, marker='^', label='Slab 3 Minimum')

# Customizing the first subplot
ax1.set_title("Plotting Agar Mins at every 15 days")
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Minimum agar values (g/l)")
ax1.set_xticks(t)  # Set x-ticks to match the time values
ax1.legend()
ax1.grid()

# Second subplot for the dimensionless time plot
ax2.plot(dim_less_time, slab_1_mins, marker='o', label='Slab 1 Minimum')
ax2.plot(dim_less_time[:len(slab_2_mins)], slab_2_mins, marker='s', label='Slab 2 Minimum')
ax2.plot(dim_less_time[:len(slab_3_mins)], slab_3_mins, marker='^', label='Slab 3 Minimum')

# Customizing the second subplot
ax2.set_title("Plotting Agar Mins (Dimensionless Time)")
ax2.set_xlabel("Dimensionless Time")
ax2.set_ylabel("Minimum agar values (g/l)")
ax2.tick_params(axis='x', rotation=45)
ax2.set_xticks(dim_less_time)  # Set x-ticks to match the time values
ax2.legend()
ax2.grid()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()


# c) Explain
print()
print("PROBLEM 2c OUTPUT")
print("i) Slab 3 takes a lot longer to get to 20 g/l compared to the other two slabs.")
print("ii) Similar to the dimensional plot, only much lower dimless time values.")
print("iii) Similarly, slab 3 takes a lot longer to get to 20 g/l compared to the other two slabs, as expected. Dimless time shouldn't change this.")
print("iv) This is a good plot showing that the aspect ratio matters greatly when we are being asked how the agar grows to 20 g/l.")

# d) Figure 2.b

# for this problem I did the exact same thing in a and b except I changed the k_m value:
# I also just got rid of the print statements below
tissue_slab_1 = [.4, .2, .01]
tissue_slab_2 = [.48, .12, .01]
tissue_slab_3 = [.56, .14, .01]

def get_mass_transfer_data_new_km(total_x, total_y):

    # defining variables
    C_e = 30 # (g/l)
    C_f = 20 # (g/l)
    D = 0.252*(10**(-4)) # (m^2/day)
    k_m = 31.536*10**(-3) # (m/day)
    dx = 0.01 # (m)
    dy = 0.01 # (m)
    
    Sh = k_m*dx/D
    interior_dt = dx**2/(4*D)
    side_dt = dx**2/(2*D*(2 + Sh))
    corner_dt = dx**2/(4*D*(1 + Sh))
    dt_maxes = [interior_dt, side_dt, corner_dt] 
    dt_max = min(dt_maxes)
    dt = round(dt_max, 3)
    Fo_m = D*dt/(dx**2)
    
    maximum = 230
    save_step = 15
    interval_steps = int(save_step / dt)
    S = int(maximum / save_step) + 1
    
    xn = int(total_x/dx) + 1
    yn = int(total_y/dy) + 1
    
    c = np.zeros((2,xn,yn)) # sets the initial conditions
    data = np.zeros((S, xn, yn))
    data[0] = c[1]
    
    # printing the first time
    # print(f"Day 0: Minimum agar = {np.min(c[1]):.4f} g/L")
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
            #print(f"Day {int(time)}: Minimum agar = {minimum:.4f} g/L")
            step_saved += 1
        
        # Check if all agar has at least target concentration
        if np.all(c[1,:,:] >= C_f):
            #print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
            break
    
    information = [data, c[1][:,1][xn-1], c[1][:,1][0], str(int(total_x*100))+"x"+str(int(total_y*100))+"x1", time, Sh]
    
    return information

slab_1_data = get_mass_transfer_data_new_km(tissue_slab_1[0], tissue_slab_1[1])
# print("x[m]: " + str(round(slab_1_data[1], 3)) + ", y[m]: " + str(round(slab_1_data[2], 3)) + ", aspect ratio: " + slab_1_data[3] + ", final time: " + str(round(slab_1_data[4],3)) + ", Sh: " + str(round(slab_1_data[5],3)))
# print()

slab_2_data = get_mass_transfer_data_new_km(tissue_slab_2[0], tissue_slab_2[1])
# print("x[m]: " + str(round(slab_2_data[1], 3)) + ", y[m]: " + str(round(slab_2_data[2], 3)) + ", aspect ratio: " + slab_2_data[3] + ", final time: " + str(round(slab_2_data[4],3)) + ", Sh: " + str(round(slab_2_data[5],3)))
# print()

slab_3_data = get_mass_transfer_data_new_km(tissue_slab_3[0], tissue_slab_3[1])
# print("x[m]: " + str(round(slab_3_data[1], 3)) + ", y[m]: " + str(round(slab_3_data[2], 3)) + ", aspect ratio: " + slab_3_data[3] + ", final time: " + str(round(slab_3_data[4],3)) + ", Sh: " + str(round(slab_3_data[5],3)))
# print()


# b) Figure 2.a
# getting the minimums of each slab at each 15 days
slab_1_mins = []
slab_2_mins = []
slab_3_mins = []

# defining the length of each array
slab_1_data_len = len(slab_1_data[0])
slab_2_data_len = len(slab_2_data[0])
slab_3_data_len = len(slab_3_data[0])

max_time_of_all_slabs = max([slab_1_data_len, slab_2_data_len, slab_3_data_len])

for i in range(slab_1_data_len-1):
    
    if np.min(slab_1_data[0][i]) == 0 and i > 0:
        break
    else:
        slab_1_mins.append(np.min(slab_1_data[0][i]))
for i in range(slab_2_data_len-1):
    if np.min(slab_2_data[0][i]) == 0 and i > 0:
        break
    else:
        slab_2_mins.append(np.min(slab_2_data[0][i]))
for i in range(slab_3_data_len-1):
    if np.min(slab_3_data[0][i]) == 0 and i > 0:
        break
    else:
        slab_3_mins.append(np.min(slab_3_data[0][i]))


t = [0, 14, 29, 44, 59, 74, 89, 104, 119, 134, 149, 164, 179, 194]

# Assuming slab_1_mins, slab_2_mins, slab_3_mins, and t are already defined

# Converting to dimensionless time
D = 0.252 * (10**(-4))  # (m^2/day)
L1 = 0.40
L2 = 0.48
L3 = 0.56

# Characteristic length
L = (L1 + L2 + L3) / 3

dim_less_time = [round(i * (D / L**2), 3) for i in t]

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# First subplot for the original time plot
ax1.plot(t, slab_1_mins, marker='o', label='Slab 1 Minimum')
ax1.plot(t[:len(slab_2_mins)], slab_2_mins, marker='s', label='Slab 2 Minimum')
ax1.plot(t[:len(slab_3_mins)], slab_3_mins, marker='^', label='Slab 3 Minimum')

# Customizing the first subplot
ax1.set_title("Plotting Agar Mins at every 15 days (New km)")
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Minimum agar values (g/l)")
ax1.set_xticks(t)  # Set x-ticks to match the time values
ax1.legend()
ax1.grid()

# Second subplot for the dimensionless time plot
ax2.plot(dim_less_time, slab_1_mins, marker='o', label='Slab 1 Minimum')
ax2.plot(dim_less_time[:len(slab_2_mins)], slab_2_mins, marker='s', label='Slab 2 Minimum')
ax2.plot(dim_less_time[:len(slab_3_mins)], slab_3_mins, marker='^', label='Slab 3 Minimum')

# Customizing the second subplot
ax2.set_title("Plotting Agar Mins (Dimensionless Time with new km)")
ax2.set_xlabel("Dimensionless Time")
ax2.set_ylabel("Minimum agar values (g/l)")
ax2.tick_params(axis='x', rotation=45)
ax2.set_xticks(dim_less_time)  # Set x-ticks to match the time values
ax2.legend()
ax2.grid()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()



# e) Explain
print("i) The agar reaches the 20 g/l cap a lot faster in the second simulation (Figure 2b), than in Figure 2a, for each aspect ratio.")
print("ii) In the dimensionless case, the Figure a and Figure b are nearly identitcal.")
print("iii) Aspect ratios greater than the critical aspect ratio may lead to instability (i.e. we will observe unrealistic behavior).")


#############
# PROBLEM 3 #
#############

# done separately