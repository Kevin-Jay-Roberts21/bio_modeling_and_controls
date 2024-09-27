#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 14:25:18 2024

@author: kevinjayroberts
"""

import matplotlib.pyplot as plt
import numpy as np

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
dT_max = (dX**2)/(2*alpha)
dT = round(dT_max, 3) 
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
plt.title('Temperature')
plt.xlabel('Dimensionless Space')
plt.ylabel('Temp in C')



# c) Output simulated U(X, T)
# plotting

# define a new u using the initial conditions
u = np.zeros((2, m))
u[0,:] = u_initial # setting the initial condition
u[1,:] = u[0,:]
loop_counter = 1

while (max(u[1,:]) >= 20):
    
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
print("PROBLEM 4 OUTPUT")
print("The final value of t is: " + str(final_time))
print()



#############
# Problem 2 #
#############

# a) Calculate U(X, T)

# b) Console ouput

# c) Describe results of decreasing Dx

# d) Plot % difference between final T and T_target (20 C) vs. the number of iterations
#    for dX = [0.1, 0.05, 0.025, 0.01]. Completely annotate the graph (labels, title, 
#    and legend with units).




#############
# Problem 3 #
#############

# a) Model heat transfer in one dimension for 1 centimeter using equation (1). Maintain
#    temp. of adjacent tissue at 37+-1C.

# b) Calculate U(X, T)

# c) Evaluate your results



