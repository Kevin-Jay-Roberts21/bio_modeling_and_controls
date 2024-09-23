#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:27:40 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt

#############
# PROBLEM 1 #
#############

# a) Simulate u(x,t) 

# initial conditions: u = 2x 0<=x<=1/2, u = 2(1-x) 1/2<=x<=1
# boundary conditions: u_x=0 = u_x=1 = 0

# defining some variables
dt = 1/1000
dx = 1/10
alpha = 1
L = 1

# space and time intervals
x_0 = 0
x_f = 1
t_0 = 0
t_f = 0.04

# creating the empty disctretized solution
n = int(t_f/dt) + 1 # iternation number
m = int(x_f/dx) + 1
u = np.zeros((n, m)) # define array 

# populating u with the initial and boundary conditions

for i in range(n):
    
    # setting the initial conditions
    if i == 0:
        for j in range(m):
            
            # setting the boundary conditions
            if j == 0:
                u[i, j] = 0
            elif j == (m-1):
                u[i, j] = 0
            else:
                # setting the boundary conditions
                if j*dx <= 0.5:     
                    u[i, j] = 2*j
                else:
                    # 2*(1 - x) changes due to the change in dx
                    u[i, j] = 2*((m-1)-j)
    else:
        for j in range(m):
            
            # setting the boundary conditions
            if j == 0:
                u[i, j] = 0
            elif j == (m-1):
                u[i, j] = 0
            else:
                u[i, j] = u[i-1,j] + dt/(dx**2)*(u[i-1,j-1] - 2*u[i-1,j] + u[i-1,j+1])

# plotting
plt.figure()
x = np.arange(0, m*dx, dx)

print(u[0,:])

t1 = int(0.01/dt)
t2 = int(0.02/dt)
t3 = int(0.03/dt)
t4 = int(0.04/dt)


plt.plot(x, u[0,:], label='Temp at t=0')
plt.plot(x, u[t1,:], label='Temp at t=0.01')
plt.plot(x, u[t2,:], label='Temp at t=0.02')
plt.plot(x, u[t3,:], label='Temp at t=0.03')
plt.plot(x, u[t4,:], label='Temp at t=0.04')


# Adding a title and labels
plt.legend()
plt.title('Initial Temp')
plt.xlabel('Space')
plt.ylabel('Temp')                











# b) Output the simulation

















#############
# PROBLEM 2 #
#############

# a) Calculate u(x,t)









# b) Output calculated u(x,t)





# c) Describe comparing simulation vs. exact solution















#############
# PROBLEM 3 #
#############

# a) Compare simulated u(x,t)















#############
# PROBLEM 4 #
#############

# a) Use equation 1 with initial and boundary conditions in Problem 1 to simulate
# u(x,t) for 0<=x<=1 at times for which u >= 0.55. Use dt_max. Inside a while loop,
# use a two-row array to compute successive values of u(x,t)










# b) Output simulated u(x,t)





# c) Console output





# d) Describe loop accuracy

