#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 13:12:33 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt


#############
# PROBLEM 1 #
#############
print("PROBELM 1 RESULTS")
# a) importing the data (these are in the same folder as Roberts_Kevin_Assign_3.py)
bacteria_data = np.loadtxt('Bacterial_Growth.csv', delimiter=',')

# b) Finding the average number of bacteria for each time step and putting data into an array
num_of_rows = len(bacteria_data)
b_growth = np.zeros((num_of_rows, 2))

for i in range(num_of_rows):
   avg = (bacteria_data[i, 1] + bacteria_data[i, 2] + bacteria_data[i, 3] + bacteria_data[i, 4])/4
   b_growth[i][1] = avg
   b_growth[i][0] = i

# c) Saving data, or average rather, to a new file called Avg_growth.txt
np.savetxt('Avg_growth.csv', b_growth, delimiter=", ")

# d) Plotting all the data and the average vs time
x = np.arange(0, num_of_rows, 1)

plt.figure()
plt.plot(x, bacteria_data[:,1], label='Bacteria Col 1')
plt.plot(x, bacteria_data[:,2], label='Bacteria Col 2')
plt.plot(x, bacteria_data[:,3], label='Bacteria Col 3')
plt.plot(x, bacteria_data[:,4], label='Bacteria Col 4')
plt.plot(x, [row[1] for row in b_growth], label='Avg. Bacteria')


# Adding a title and labels
plt.title('Plotting all Bacterias and the Avg. Growth')
plt.xlabel('Time in Hours')
plt.ylabel('Bacteria Count')

# e) Outputting all of the growth and average to the consol
np.set_printoptions(precision=3, threshold=100, edgeitems=3, suppress=False)
all_data = np.zeros((num_of_rows, 6))

for i in range(num_of_rows):
    for j in range(5):
        all_data[i][j] = bacteria_data[i][j] 
    all_data[i][5] = b_growth[i][1]

print("Here is the outputted data from everything, includeing time, all bacterias and average")
print("time (hrs), bacteria 1, bacteria 2, bacteria 3, bacteria 4, average")
print(all_data)
print()
print()

#############
# PROBLEM 2 #
#############
print("PROBLEM 2 RESULTS")

# a) Use a for loop to compute value of S, I and R at each time step in one array

S_0 = 666666
I_0 = 1
R_0 = 0
N = 666667 # total population
final_t = 90
dt = 0.5
n = int(final_t/dt) + 1
b = 3/5
k = 1/6

total_pop = np.zeros((n, 5)) # includes i, time, S, I and R respectively

# initialization
total_pop[0, 2] = S_0
total_pop[0, 3] = I_0
total_pop[0, 4] = R_0
time = 0

for i in range(1, n):
    total_pop[i, 0] = i
    time += dt
    total_pop[i, 1] += time
    total_pop[i, 2] = total_pop[i-1, 2] - b/N*dt*total_pop[i-1, 2]*total_pop[i-1, 3] # updating S
    total_pop[i, 3] = total_pop[i-1, 3] + dt*(b/N*total_pop[i-1, 2]*total_pop[i-1, 3] - k*total_pop[i-1, 3]) # updating I    
    total_pop[i, 4] = total_pop[i-1, 4] + k*dt*total_pop[i-1, 3] # updating R

# b) Outputting plot data and computed array to the console

x = np.arange(0, n*dt, dt)

plt.figure()
plt.plot(x, [row[2] for row in total_pop], marker='o', label='Population S')
plt.plot(x, [row[3] for row in total_pop], marker='o', label='Population I')
plt.plot(x, [row[4] for row in total_pop], marker='o', label='Population R')

# Adding a title and labels
plt.title('Plotting Susceptibles, Infecteds and Recoverd')
plt.xlabel('Time in Days')
plt.ylabel('Population')

np.set_printoptions(precision=3, threshold=100, edgeitems=3, suppress=False)
print("Data from the SIR Model")
print("i, time (days), Susceptibles, Infected, Recovered")
print(total_pop)
print()
print()

#############
# PROBLEM 3 #
#############
print("PROBLEM 3 RESULTS")

# a) Develope a finite difference simulation for the given equations

h = 10 # effect of presence of prey on predator
k = 10000 # carrying capacity of hares
r = 0.9 # per capita rate of prey increase
S = 0.3 # per capita rate of predator increase
B = 1000 # normalization factor
dt = 0.25 # time increment in years
total_time = 100 # in years
n = int(total_time/dt) + 1
N_0 = 1000
P_0 = 100

# defining a function that returns the predator prey data given a value of A
def get_pred_prey_data(A):

    pred_prey_array = np.zeros((n, 4)) # iteration, time, prey, predator
    pred_prey_array[0, 2] = N_0
    pred_prey_array[0, 3] = P_0
    time = 0
    
    for i in range(1, n):
        pred_prey_array[i, 0] = i
        time += dt
        pred_prey_array[i, 1] += time
        pred_prey_array[i, 2] = pred_prey_array[i-1, 2] + dt*pred_prey_array[i-1, 2]*(r*(1 - pred_prey_array[i-1, 2]/k) - A*pred_prey_array[i-1, 2]*pred_prey_array[i-1, 3]/(pred_prey_array[i-1, 2]**2 + B**2)) # updating Prey
        pred_prey_array[i, 3] = pred_prey_array[i-1, 3] + dt*pred_prey_array[i-1, 3]*(S*(1 - h*pred_prey_array[i-1, 3]/pred_prey_array[i-1, 2])) # updating Predator    

    return pred_prey_array

pred_prey_A_is_8 = get_pred_prey_data(8) 


# b) plotting the data
x = np.arange(0, n*dt, dt)

plt.figure()

fig, ((ax1,ax2)) = plt.subplots(1,2)
fig.suptitle("Plotting Predator and Prey Populations")
ax1.plot(x, [row[2] for row in pred_prey_A_is_8], c="blue", marker="D", label="Prey Population") #
ax2.plot(x, [row[3] for row in pred_prey_A_is_8], c="red", marker="^", label="Predator Population")
ax1.set_xlabel('Time in Years')
ax1.set_ylabel('Prey Count')
ax2.set_xlabel('Time in Years')
ax2.set_ylabel('Predator Count')
ax1.set_title("Prey Dynamics") # Plot title
ax2.set_title("Predator Dynamics")
plt.tight_layout()
plt.show()

# c) Expand code to compute and plot N(t) and P(t) for the 5 different A values

# getting all the data for each value of A
pred_prey_A_is_2 = get_pred_prey_data(2)
pred_prey_A_is_4 = get_pred_prey_data(4)
pred_prey_A_is_5 = get_pred_prey_data(5)
pred_prey_A_is_6 = get_pred_prey_data(6)

# creating the figure
x = np.arange(0, n*dt, dt)

plt.figure()

fig, ((ax1,ax2)) = plt.subplots(1,2)
fig.suptitle("Plotting Predator and Prey Populations")
ax1.plot(x, [row[2] for row in pred_prey_A_is_2], c="blue", marker="D", label="Prey Population A = 2")
ax1.plot(x, [row[2] for row in pred_prey_A_is_4], c="green", marker="^", label="Prey Population A = 4")
ax1.plot(x, [row[2] for row in pred_prey_A_is_5], c="orange", marker="o", label="Prey Population A = 5")
ax1.plot(x, [row[2] for row in pred_prey_A_is_6], c="purple", marker="s", label="Prey Population A = 6")
ax1.plot(x, [row[2] for row in pred_prey_A_is_8], c="red", marker="*", label="Prey Population A = 8")
ax2.plot(x, [row[3] for row in pred_prey_A_is_2], c="blue", marker="D", label="Predator Population A = 2")
ax2.plot(x, [row[3] for row in pred_prey_A_is_4], c="green", marker="^", label="Predator Population A = 4")
ax2.plot(x, [row[3] for row in pred_prey_A_is_5], c="orange", marker="o", label="Predator Population A = 5")
ax2.plot(x, [row[3] for row in pred_prey_A_is_6], c="purple", marker="s", label="Predator Population A = 6")
ax2.plot(x, [row[3] for row in pred_prey_A_is_8], c="red", marker="*", label="Predator Population A = 8")
ax1.set_xlabel('Time in Years')
ax1.set_ylabel('Prey Count')
ax2.set_xlabel('Time in Years')
ax2.set_ylabel('Predator Count')
ax1.set_title("Prey Dynamics") # Plot title
ax2.set_title("Predator Dynamics")
plt.tight_layout()
plt.show()



