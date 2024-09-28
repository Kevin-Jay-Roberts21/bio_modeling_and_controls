#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 14:25:18 2024

@author: kevinjayroberts
"""

import matplotlib.pyplot as plt
import numpy as np

def perc_diff(A, B):
    return np.abs(A-B)/((A+B)/2)*100

# setting printing params
np.set_printoptions(precision=3, edgeitems=3, suppress=True)

#############
# Problem 1 #
#############
print()
print("PROBLEM 1a OUTPUT")

# a) Parameters

dX = 0.1 # (meter)
k = 54 # (W/m-C)
p = 780 # (kg/m^3)
Cp = 460 # (J/kg-C)
alpha = k/(p*Cp)
dT_max = (dX**2)/(2*alpha) # in seconds
dT = round(dT_max, 3) # in seconds
print("The following values for alpha, dT_max and dT are: ")
print("alpha = " + str(round(alpha, 5)))
print("dT_max = " + str(round(dT_max, 5)))
print("dT = " + str(dT))
print("The since alpha is much greater than 1 in this case, this means that the heat")
print("will diffuse much faster than what it did in the last problem. Additionally,")
print("the dT in this example will increase since our criteria for dT_max changes")
print("due to our alpha value.")

# b) Initial and Boundary Conditions

# plotting the initial conditions
x_f = 1
m = int(x_f/dX) + 1
u_initial = np.zeros((1, m))
for j in range(m):
    # setting the boundary conditions
    if j == 0:
        u_initial[0, j] = 0
    elif j == (m-1):
        u_initial[0, j] = 0
    else:
        u_initial[0, j] = 100

plt.figure()
x = np.arange(0, m*dX, dX)

plt.plot(x, u_initial[0,:], 'o', label='Temp at t=0')

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Temperature Diffusion of Rod')
plt.xlabel('Meters')
plt.ylabel('Temp in C')



# c) Output simulated U(X, T)
# plotting

# define a new u using the initial conditions
u = np.zeros((2, m))
u[0,:] = u_initial # setting the initial condition
u[1,:] = u[0,:]
loop_counter = 1

mid = int(0.5/dX)

while (u[1, mid] >= 20):
    
    u[0,:] = u[1,:]
    
    for j in range(m):
        
        # setting the boundary conditions
        if j == 0:
            u[1, j] = 0
        elif j == (m-1):
            u[1, j] = 0
        else:
            u[1, j] = u[0,j] + dT/(dX**2)*alpha*(u[0,j-1] - 2*u[0,j] + u[0,j+1])
    
    
    # plotting only the every 10 iteration of the while loop
    if loop_counter % 10 == 0:
        
        plt.plot(x, u[1,:], 'o', label='t={}, u_max = {}'.format(round(loop_counter*dT,3), round(max(u[1,:]),3)))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    loop_counter += 1

final_time = (loop_counter-1)*dT
plt.text(1.05, 0.75, 'The final time is ' + str(final_time) + " secs", transform=plt.gca().transAxes, fontsize=10, color='black')
print("Percetange difference between final Temp and 20 C at the middle of the rod is: " + str(round(perc_diff(u[1,mid], 20), 3)))
print("The number of iterations to reach this U_max was: " + str(loop_counter-1))












#############
# Problem 2 #
#############

# a) Calculate U(X, T)

print()
print("PROBLEM 2a OUTPUT")

dX = 0.05 # (meter)
k = 54 # (W/m-C)
p = 780 # (kg/m^3)
Cp = 460 # (J/kg-C)
alpha = k/(p*Cp)
dT_max = (dX**2)/(2*alpha) # in seconds
dT = round(dT_max, 3) # in seconds
print("The following values for alpha, dT_max and dT are: ")
print("alpha = " + str(round(alpha, 5)))
print("dT_max = " + str(round(dT_max, 5)))
print("dT = " + str(dT))

# plotting the initial conditions
x_f = 1
m = int(x_f/dX) + 1
u_initial = np.zeros((1, m))
for j in range(m):
    # setting the boundary conditions
    if j == 0:
        u_initial[0, j] = 0
    elif j == (m-1):
        u_initial[0, j] = 0
    else:
        u_initial[0, j] = 100

plt.figure()
x = np.arange(0, m*dX, dX)

plt.plot(x, u_initial[0,:], 'o', label='Temp at t=0')

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Temperature Diffusion of Rod')
plt.xlabel('Meters')
plt.ylabel('Temp in C')


# define a new u using the initial conditions
u = np.zeros((2, m))
u[0,:] = u_initial # setting the initial condition
u[1,:] = u[0,:]
loop_counter = 1

print(u)

mid = int(0.5/dX)

while (u[1, mid] >= 20):
    
    u[0,:] = u[1,:]
    
    for j in range(m):
        
        # setting the boundary conditions
        if j == 0:
            u[1, j] = 0
        elif j == (m-1):
            u[1, j] = 0
        else:
            u[1, j] = u[0,j] + dT/(dX**2)*alpha*(u[0,j-1] - 2*u[0,j] + u[0,j+1])
    
    
    # plotting only the every 10 iteration of the while loop
    if loop_counter % 40 == 0:
        
        plt.plot(x, u[1,:], 'o', label='t={}, u_max = {}'.format(round(loop_counter*dT,3), round(max(u[1,:]),3)))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    loop_counter += 1

final_time = round((loop_counter-1)*dT, 3)
plt.text(1.05, 0.75, 'The final time is ' + str(final_time) + " secs", transform=plt.gca().transAxes, fontsize=10, color='black')
print("Percetange difference between final Temp and 20 C at the middle of the rod is: " + str(round(perc_diff(u[1,mid], 20), 3)))
print()
# b) Console ouput

print("PROBLEM 2b OUTPUT")
print("The final value of t is: " + str(final_time) + " seconds.")
print("The final temperature across the rod is at this time is: ")
print(u[1,:])
print("The number of iterations to reach this U_max was: " + str(loop_counter-1))
print()

# c) Describe results of decreasing dX
print("PROBLEM 2c OUTPUT")
print("i) The dT_max decreased significantly as dX decreased. And the temperature")
print("decayed 'smoother' as dX decreased.\n")
print("ii) As dX decreased, the percentage difference also decreased which is good news.")
print("This means that we gave an approximation that got closer to what we expect.")
print("The value of dT also dexreased.\n")
print("iii) The number of loop iterations increased significantly. This makes\n sense since we decreased dX.")
print()




# d) Plot % difference between final T and T_target (20 C) vs. the number of iterations
#    for dX = [0.1, 0.05, 0.025, 0.01]. Completely annotate the graph (labels, title, 
#    and legend with units).

# defining a function to find u given a dX value
def approx_u(dX):
    
    k = 54 # (W/m-C)
    p = 780 # (kg/m^3)
    Cp = 460 # (J/kg-C)
    alpha = k/(p*Cp)
    dT_max = (dX**2)/(2*alpha) # in seconds
    dT = round(dT_max, 3) # in seconds
    x_f = 1
    m = int(x_f/dX) + 1
    u_initial = np.zeros((1, m))
    
    for j in range(m):
        # setting the boundary conditions
        if j == 0:
            u_initial[0, j] = 0
        elif j == (m-1):
            u_initial[0, j] = 0
        else:
            u_initial[0, j] = 100

    # define a new u using the initial conditions
    u = np.zeros((2, m))
    u[0,:] = u_initial # setting the initial condition
    u[1,:] = u[0,:]
    loop_counter = 1

    mid = int(0.5/dX)

    while (u[1, mid] >= 20):
        
        u[0,:] = u[1,:]
        
        for j in range(m):
            
            # setting the boundary conditions
            if j == 0:
                u[1, j] = 0
            elif j == (m-1):
                u[1, j] = 0
            else:
                u[1, j] = u[0,j] + dT/(dX**2)*alpha*(u[0,j-1] - 2*u[0,j] + u[0,j+1])
        
        loop_counter += 1

    u_info = [loop_counter-1, perc_diff(u[1,mid], 20)] # returning percentage difference and number of iterations

    return u_info

dX_values = [0.1, 0.05, 0.025, 0.01]
u_dX_0 = approx_u(dX_values[0])
u_dX_1 = approx_u(dX_values[1])
u_dX_2 = approx_u(dX_values[2])
u_dX_3 = approx_u(dX_values[3])

# plotting
plt.figure()
plt.xlim([0, 4000])
plt.ylim([0, 5])

x = [u_dX_0[0], u_dX_1[0], u_dX_2[0], u_dX_3[0]]
y = [u_dX_0[1], u_dX_1[1], u_dX_2[1], u_dX_3[1]]
labels = ["dX=0.1", "dX=0.05", "dX=0.025", "dX=0.01"]

plt.scatter(x, y, color='blue')

for i in range(len(x)):
    plt.text(x[i]+ 50, y[i], labels[i], fontsize=12, ha='left', va='bottom')

# Adding a title and labels
plt.title('Analyzing Percentage Difference with different dX\'s')
plt.xlabel('Number of iterations')
plt.ylabel('Percentage Difference')






#############
# Problem 3 #
#############
print()


# a) Model heat transfer in one dimension for 1 centimeter using equation (1). Maintain
#    temp. of adjacent tissue at 37+-1C.

# b) Calculate U(X, T)

# c) Evaluate your results



