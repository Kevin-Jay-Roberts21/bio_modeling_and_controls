#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:27:40 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt

# setting printing params
np.set_printoptions(precision=3, edgeitems=3, suppress=True)

#############
# PROBLEM 1 #
#############

# a) Simulate u(x,t) 

# initial conditions: u = 2x 0<=x<=1/2, u = 2(1-x) 1/2<=x<=1
# boundary conditions: u_x=0 = u_x=1 = 0

# populating u with the initial and boundary conditions

def approx(dx, dt):
    
    # space and time intervals
    x_0 = 0
    x_f = 1
    t_0 = 0
    t_f = 0.04
    
    # creating the empty disctretized solution
    n = int(t_f/dt) + 1 # iternation number
    m = int(x_f/dx) + 1
    u = np.zeros((n, m)) # define array 
    
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
                        u[i, j] = 2*j*dx
                    else:
                        # 2*(1 - x) changes due to the change in dx
                        u[i, j] = 2*((m-1)-j)*dx
        else:
            for j in range(m):
                
                # setting the boundary conditions
                if j == 0:
                    u[i, j] = 0
                elif j == (m-1):
                    u[i, j] = 0
                else:
                    u[i, j] = u[i-1,j] + dt/(dx**2)*(u[i-1,j-1] - 2*u[i-1,j] + u[i-1,j+1])
    return u

# defining some variables
dt = 1/1000
dx = 1/10
u = approx(dx, dt)

# b) Output the simulation
# plotting
plt.figure()
x = np.arange(0, len(u[0,:])*dx, dx)

t1 = int(0.01/dt)
t2 = int(0.02/dt)
t3 = int(0.03/dt)
t4 = int(0.04/dt)

plt.plot(x, u[0,:], label='Temp at t=0, u_max = {}'.format(round(max(u[0,:]),3)))
plt.plot(x, u[t1,:], label='Temp at t=0.01, u_max = {}'.format(round(max(u[t1,:]),3)))
plt.plot(x, u[t2,:], label='Temp at t=0.02, u_max = {}'.format(round(max(u[t2,:]),3)))
plt.plot(x, u[t3,:], label='Temp at t=0.03, u_max = {}'.format(round(max(u[t3,:]),3)))
plt.plot(x, u[t4,:], label='Temp at t=0.04, u_max = {}'.format(round(max(u[t4,:]),3)))

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Approx. Temperature Diffusion Profile')
plt.xlabel('Dimensionless Space')
plt.ylabel('Dimensionless Temp')


#############
# PROBLEM 2 #
#############

# a) Calculate u(x,t)

# finding the exact solution

def exact(dx, dt):
    
    # space and time intervals
    x_0 = 0
    x_f = 1
    t_0 = 0
    t_f = 0.04
    
    # creating the empty disctretized solution
    n = int(t_f/dt) + 1 # iternation number
    m = int(x_f/dx) + 1
    e = np.zeros((n, m)) # define array 
    
    N = 100
    for i in range(n):
        
        # setting the initial conditions
        if i == 0:
            for j in range(m):
                
                # setting the boundary conditions
                if j == 0:
                    e[i, j] = 0
                elif j == (m-1):
                    e[i, j] = 0
                else:
                    # setting the boundary conditions
                    if j*dx <= 0.5:     
                        e[i, j] = 2*j*dx
                    else:
                        # 2*(1 - x) changes due to the change in dx
                        e[i, j] = 2*((m-1)-j)*dx
        else:
            for j in range(m):
                
                # setting the boundary conditions
                if j == 0:
                    e[i, j] = 0
                elif j == (m-1):
                    e[i, j] = 0
                else:
                    exact_sum = 0
                    for k in range(1, N+1):
                        exact_sum += 1/(k**2)*np.sin(k*np.pi/2)*np.sin(n*np.pi*j*dx)*np.exp(-(k**2)*(np.pi**2)*i*dt)
                    e[i, j] = 8/(np.pi**2)*exact_sum
    
    return e
    

# defining some variables
dt = 1/1000
dx = 1/10
e = exact(dx, dt) # getting the exact solution

# b) Output calculated u(x,t)
# plotting
plt.figure()
x = np.arange(0, len(e[0,:])*dx, dx)

t1 = int(0.01/dt)
t2 = int(0.02/dt)
t3 = int(0.03/dt)
t4 = int(0.04/dt)

plt.plot(x, e[0,:], label='Temp at t=0, u_max = {}'.format(round(max(e[0,:]),3)))
plt.plot(x, e[t1,:], label='Temp at t=0.01, u_max = {}'.format(round(max(e[t1,:]),3)))
plt.plot(x, e[t2,:], label='Temp at t=0.02, u_max = {}'.format(round(max(e[t2,:]),3)))
plt.plot(x, e[t3,:], label='Temp at t=0.03, u_max = {}'.format(round(max(e[t3,:]),3)))
plt.plot(x, e[t4,:], label='Temp at t=0.04, u_max = {}'.format(round(max(e[t4,:]),3)))

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Exact Temperature Diffusion Profile')
plt.xlabel('Dimensionless Space')
plt.ylabel('Dimensionless Temp')



# c) Describe comparing simulation vs. exact solution
print()
print("PROBLEM 2 OUTPUT")
print("The exact solution tends to diffuse faster and 'smoother' than the approximation.")
print("To get the approximated solution closer to the exact solution we can decrease the")
print("time and space steps dt and dx.")
print()
print()



#############
# PROBLEM 3 #
#############

# a) Compare simulated u(x,t)
dx = 0.1
dt1 = (dx**2)/2 # dt right at the threshold of stability
dt2 = 1.1*((dx**2)/2) # unstable dt

u_dt1 = approx(dx, dt1)
u_dt2 = approx(dx, dt2)

# plotting
# plotting
plt.figure()
x = np.arange(0, len(u_dt1[0,:])*dx, dx)

t1 = int(0.01/dt1)
t2 = int(0.02/dt1)
t3 = int(0.03/dt1)
t4 = int(0.04/dt1)

plt.plot(x, u_dt1[0,:], label='t=0, dt=0.005, u_max = {}'.format(round(max(u_dt1[0,:]),3)))
plt.plot(x, u_dt1[t1,:], label='t=0.01, dt=0.005, u_max = {}'.format(round(max(u_dt1[t1,:]),3)))
plt.plot(x, u_dt1[t2,:], label='t=0.02, dt=0.005, u_max = {}'.format(round(max(u_dt1[t2,:]),3)))
plt.plot(x, u_dt1[t3,:], label='t=0.03, dt=0.005, u_max = {}'.format(round(max(u_dt1[t3,:]),3)))
plt.plot(x, u_dt1[t4,:], label='t=0.04, dt=0.005, u_max = {}'.format(round(max(u_dt1[t4,:]),3)))
plt.plot(x, u_dt2[0,:], label='t=0, dt=0.0055, u_max = {}'.format(round(max(u_dt2[0,:]),3)))
plt.plot(x, u_dt2[t1,:], label='t=0.01, dt=0.0055, u_max = {}'.format(round(max(u_dt2[t1,:]),3)))
plt.plot(x, u_dt2[t2,:], label='t=0.02, dt=0.0055, u_max = {}'.format(round(max(u_dt2[t2,:]),3)))
plt.plot(x, u_dt2[t3,:], label='t=0.03, dt=0.0055, u_max = {}'.format(round(max(u_dt2[t3,:]),3)))
plt.plot(x, u_dt2[t4,:], label='t=0.04, dt=0.0055, u_max = {}'.format(round(max(u_dt2[t4,:]),3)))

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Comparison of Approx. Sol. with diff. dt values')
plt.xlabel('Dimensionless Space')
plt.ylabel('Dimensionless Temp')


# b) Describe sensitivity of dt
print("PROBLEM 3 OUTPUT")
print("When we pick a dt value larger than (dx^2)/2, we see some instability in")
print("the system, i.e. the temperature does not diffuse like we exact it to. When")
print("pick a dt value lower than this threshold, we see stability in the system.")
print()
print()

#############
# PROBLEM 4 #
#############

# a) Use equation 1 with initial and boundary conditions in Problem 1 to simulate
# u(x,t) for 0<=x<=1 at times for which u >= 0.55. Use dt_max. Inside a while loop,
# use a two-row array to compute successive values of u(x,t)

# space and time intervals
x_0 = 0
x_f = 1
t_0 = 0
t_f = 0.04

dt = 0.005
dx = 1/10

# creating the empty disctretized solution
n = int(t_f/dt) + 1 # iternation number
m = int(x_f/dx) + 1
u = np.zeros((2, m)) # define array 

# initializing the array
for j in range(m):
    
    # setting the boundary conditions
    if j == 0:
        u[0, j] = 0
    elif j == (m-1):
        u[0, j] = 0
    else:
        # setting the boundary conditions
        if j*dx <= 0.5:     
            u[0, j] = 2*j*dx
        else:
            # 2*(1 - x) changes due to the change in dx
            u[0, j] = 2*((m-1)-j)*dx

plt.figure()
x = np.arange(0, len(u[0,:])*dx, dx)

plt.plot(x, u[0,:], label='t={}, u_max = {}'.format(0, round(max(u[0,:]),3)))

t1 = int(0.01/dt1)
t2 = int(0.02/dt1)
t3 = int(0.03/dt1)
t4 = int(0.04/dt1)

u[1,:] = u[0,:]

loop_counter = 1
while (max(u[1,:]) >= 0.55):
    
    u[0,:] = u[1,:]
    
    for j in range(m):
        
        # setting the boundary conditions
        if j == 0:
            u[1, j] = 0
        elif j == (m-1):
            u[1, j] = 0
        else:
            u[1, j] = u[0,j] + dt/(dx**2)*(u[0,j-1] - 2*u[0,j] + u[0,j+1])
    
    
    
    plt.plot(x, u[1,:], label='t={}, u_max = {}'.format(round(loop_counter*dt,3), round(max(u[1,:]),3)))
    # Adding a title and labels
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Plotting Until u_max decreases to 0.55')
    plt.xlabel('Dimensionless Space')
    plt.ylabel('Dimensionless Temp')
    
    loop_counter += 1


# b) Output simulated u(x,t)
# This was done in the while loop

# c) Console output
final_time = (loop_counter-1)*dt
print("PROBLEM 4 OUTPUT")
print("The final value of t is: " + str(final_time))
print("The final value of u along the length is: ")
print(u[1,:])
print()

# d) Describe loop accuracy
print("The simulation did end at the target value of u because the max of u at")
print("the final iteration is less than 0.55. If it didn't, we could decrease dt.")


