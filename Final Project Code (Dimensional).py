# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 17:14:38 2024

@author: Kevin Roberts
"""

import matplotlib.pyplot as plt
import numpy as np

# defining all of the constant parameters

v = 2 # (dimless)
n_0 = (1/(10**-5))**3 # (m^3)
c_0 = 6.6*(10**-7) # (mol/m^3)
k = 6.31*(10**-3) # (1/h)
L = 0.002 # (m)
T = 25 # (h) # they run their experiment for 37 hours, we're doing 25 hours.
D_c = 25*k*(L**2) # (m^2/h)
mu = 5.75*(10**-20) # (mol*h/cell)
c_hat = 2*(10**-6) # (mol/m^3)
A = 560 # (dimless)
sigma = 4000 # (dimless)
delta = np.log(2) # (1/h)

alpha_star = 0.01 # (dimless)
beta_star = 0.1 # (dimless)
alpha = alpha_star*k*(L**2) # (m^2/h)
beta = beta_star*k*(L**2) # (m^2/h)

alpha_1_star = 0.9 # (dimless)
beta_1_star = 0.1 # (dimless)
alpha_1 = 0.9*k # (1/h)
beta_1 = 0.1*k # (1/h)

# defining the functions:
def D_n(c):
    return alpha*c/c_0 + beta

def s(c):
    return alpha_1*c/c_0 + beta_1

def B(n):
    if n < 0.2:
        return sigma
    elif 0.2 <= n <= 0.4:
        return sigma*(2 - 5*n)
    else: # if n/n_0 > 0.4
        return 0

def f(n):
    return 2*(k*c_0*A + k*c_0*B(n)) # Note that the 2 multiplier isn't in the original equation. It doubles cell production. Manually changed for visuals.

def percent_difference(val1, val2):

    return round(2 * (val1 - val2) / (val1 + val2), 3) * 100

def Model(c_0):
    # defining the time and space step as well as total length and total time
    dx = 0.0001 # (m)
    dt = 0.001 # (hrs) checking the status every hour
    total_length = L
    total_time = T
    n = int(total_time/dt) + 1 # number of time steps
    m = int(2*total_length/dx) + 1 # number of space steps

    # defining the cell density and EGF concentration arrays
    n_data = np.zeros((n, m))
    c_data = np.zeros((n, m))

    n_snapshot = []
    c_snapshot = []
    time_list = []

    # Running the loop
    # time step

    num_sections = 9  # Change this to the number of sections you want
    t_hours = np.linspace(1, 25, num_sections).tolist()  # Generates evenly spaced values
    hours = [int(1000 * i) for i in t_hours] # The time steps that will be saved. 
    
    # Time step
    for i in range(n):
        # distance steps
        for j in range(m):
            
            if i == 0:
                
                # initial conditions
                if 0 <= j < (L/2)/dx: # (-inf < x < 0) 
                    n_data[i, j] = n_0
                    c_data[i, j] = c_0
                else: # (0 <= x < L)
                    n_data[i, j] = 0
                    c_data[i, j] = 0
            else:
                
                # setting the boundary condition
                if j == 0: # (x = -inf, left end) 
                    n_data[i, j] = n_0
                    c_data[i, j] = c_0
                
                elif j == m-1: # (neumann boundary conditions, right end)
                
                    # neumann boundary conditions
                    n_data[i, j] = n_data[i, j-1]
                    c_data[i, j] = c_data[i, j-1]
                
                    # # dirichlet boundary conditions
                    # n_data[i, j-1] = 0
                    # c_data[i, j-1] = 0
                else:
                    
                    # another approach (centered difference for the mitigation term)
                    n_data[i, j] = n_data[i-1, j] + dt*(1/(2*dx**2)*((D_n(c_data[i-1, j+1]) + D_n(c_data[i-1, j]))*(n_data[i-1, j+1] - n_data[i-1, j]) - (D_n(c_data[i-1, j]) + D_n(c_data[i-1, j-1]))*(n_data[i-1, j] - n_data[i-1, j-1])) + s(c_data[i-1, j])*n_data[i-1, j]*(v - n_data[i-1, j]/n_0) - k*n_data[i-1,j])
                    
                    
                    # leap frogging result (forward difference for the mitigation term)
                    # n_data[i, j] = n_data[i-1, j] + dt*(1/(2*(dx**2))*D_n(c_data[i-1, j])*(n_data[i-1, j+1] - n_data[i-1, j-1]) + s(c_data[i-1, j])*n_data[i-1, j]*(v - n_data[i-1, j]/n_0) - k*n_data[i-1,j])
                    
                    
                    c_data[i, j] = c_data[i-1, j] + dt*(D_c*(c_data[i-1, j+1] - 2*c_data[i-1, j] + c_data[i-1, j+1])/(dx**2) + f(n_data[i-1, j]) - mu*c_data[i-1, j]/(c_hat + c_data[i-1, j])*n_data[i-1, j] - delta*c_data[i-1, j])

    # Pass through the calculated concentrations (both c and n) and take a "snapshot" of the data at desired print times.
    for i in range(len(c_data)):
        if i in hours or (i == 0):
            c_snapshot.append(np.copy(c_data[i]))
            n_snapshot.append(np.copy(n_data[i]))
            time_list.append(i*dt)


    # b) Output the simulation
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))
 
    x = np.linspace(-L/2, L/2, len(c_data[0,:]))
    
    for i in range(len(c_snapshot)):
        axs[0].plot(x, n_snapshot[i], label = "Time = " + str(time_list[i]) + " hours")
        axs[1].plot(x, c_snapshot[i], label = "Time = " + str(time_list[i]) + " hours")

    title_font = 30
    axis_font = 20
    tick_size = 10

    axs[0].set_title("Cell Density Spatial Profile", fontsize = title_font)
    axs[0].set_ylabel(r"Cell Density ($\frac{cell}{m^3}$)", fontsize = axis_font)
    axs[0].set_xlabel("Wound Profile (m)", fontsize = axis_font)
    axs[0].tick_params(axis='both', which='major', labelsize=tick_size)
    axs[0].legend()

    axs[1].set_title("EGF Concentration Spatial Profile", fontsize = title_font)
    axs[1].set_ylabel(r"EGF Concentration ($\frac{mol}{m^3}$)", fontsize = axis_font)
    axs[1].set_xlabel("Wound Profile (m)", fontsize = axis_font)
    axs[1].tick_params(axis='both', which='major', labelsize=tick_size)
    axs[1].legend()

    plt.tight_layout(pad=3.0)
    plt.show()
    return min(n_data[-1])
    
print(Model(c_0))