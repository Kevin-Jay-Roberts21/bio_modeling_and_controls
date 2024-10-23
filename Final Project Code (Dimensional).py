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
c_0 = 6.6*(10**-10) # (M)
k = 6.31*(10**-3) # (1/h)
L = 0.004 # (m)
T = 37 # (h) # they run their experiment for 37 hours
D_c = 25*k*(L**2) # (m^2/h)
mu = 5.75*(10**-20) # (mol*h/cell)
c_hat = 2*(10**-9) # (M)
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
    if n/n_0 < 0.2:
        return sigma
    elif 0.2 <= n/n_0 <= 0.4:
        return sigma*(2 - 5*n/n_0)
    else: # if n/n_0 > 0.4
        return 0

def f(n):
    return k*c_0*A + k*c_0*B(n)

# defining the time and space step as well as total length and total time
dx = 0.0001 # (m)
dt = 1 # (hrs) checking the status every hour
total_length = L
total_time = T
n = int(total_time/dt) + 1 # number of time steps
m = int(total_length/dx) + 1 # number of space steps

# defining the boundary and initial conditions

# initial conditions
n_time_forpositivex_0 = 0 # initial condition (for 0 <= x < L)
c_time_forpositivex_0 = 0 # initial condition (for 0 <= x < L)
n_time_fornegativex_0 = n_0 # initial condition (for -inf < x < 0)
c_time_fornegativex_0 = c_0 # initial condition (for -inf < x < 0)

# boundary conditions
n_left_boudnary = n_0 # boundary condition (n(-inf, t) = n_0 for any t>=0)
c_left_boundary = c_0 # boundary condition (c(-inf, t) = c_0 for any t>=0)
n_right_boundary = 0 # neumann condition (n_x(L,t) = 0 for any t>=0)
c_right_boundary = 0 # neumann condition (c_x(L,t) = 0 for any t >= 0)

# defining the cell density and EGF concentration arrays
n_data = np.zeros((n, m))
c_data = np.zeros((n, m))

# Running the loop
for i in range(n):
    for j in range(m):
        
        
        if i == 0 and 0 <= j < (L/2)/dx:
            
            if j == 0: # left boundary conditions
                n_data[i, j] = n_0
                c_data[i, j] = c_0
            else: # initial conditions
                n_data[i, j] = n_0 # (-inf < x < 0)
                c_data[i, j] = c_0 # (-inf < x < 0)
        elif i == 0 and (L/2)/dx < j < L/dx:
            n_data[i, j] = 0
            c_data[i, j] = 0
        
        
        
        
        # the regular discretization
        else:
            n_data[i, j] = n_data(i-1, j) + dt*(1/(2*(dx**2))*D_n(c_data[i-1, j])*(n_data[i-1, j+1] - n_data[i-1, j-1]) + s(c_data[i-1, j])*n_data[i-1, j]*(v - n_data[i-1, j]/n_0) - k*n_data[i-1,j])
            c_data[i, j] = c_data(i-1, j) + dt*(D_c/(dx**2)*(c_data[i-1, j+1] - 2*c_data[i-1, j] + c[i-1, j+1]) + f(n_data(i-1, j)) - mu*c_data[i-1, j]/(c_hat + c_data[j-1, i])*n_data[i-1, j] - delta*c_data[i-1, j])













# defining the EGF concentration and the density
c_data = np.zeros(())