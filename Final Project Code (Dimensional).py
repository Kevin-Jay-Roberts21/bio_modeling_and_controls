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
c_0 = 6.6*2*(10**-7) # (mol/m^3)
k = 6.31*(10**-3) # (1/h)
L = 0.002 # (m)
T = 37 # (h) # they run their experiment for 37 hours
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
    return k*c_0*A + k*c_0*B(n)

# defining the time and space step as well as total length and total time
dx = 0.0001 # (m)
dt = 0.001 # (hrs) checking the status every hour
total_length = L
total_time = T
n = int(total_time/dt) + 1 # number of time steps
m = int(total_length/dx) + 1 # number of space steps

# defining the cell density and EGF concentration arrays
n_data = np.zeros((n, m))
c_data = np.zeros((n, m))

# Running the loop
for i in range(n):
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
        

# b) Output the simulation
# plotting
plt.figure()
x = np.linspace(-L/2, L/2, len(n_data[0,:]))

t1_hours = 7
t2_hours = 8
t3_hours = 9
t4_hours = 9
t5_hours = 12
t6_hours = 15
t7_hours = 18
t8_hours = 21
t9_hours = 25


t0 = 0
t1 = int(t1_hours/dt)
t2 = int(t2_hours/dt)
t3 = int(t3_hours/dt)
t4 = int(t4_hours/dt)
t5 = int(t5_hours/dt)
t6 = int(t6_hours/dt)
t7 = int(t7_hours/dt)
t8 = int(t8_hours/dt)
t9 = int(t9_hours/dt)

#plt.plot(x, n_data[t0,:], label=f'Density at t={t0} hours')
# plt.plot(x, n_data[t1,:], label=f'Density at t={t1_hours} hours')
# plt.plot(x, n_data[t2,:], label=f'Density at t={t2_hours} hours')
# plt.plot(x, n_data[t3,:], label=f'Density at t={t3_hours} hours')
plt.plot(x, n_data[t4,:], label=f'Density at t={t4_hours} hours')
plt.plot(x, n_data[t5,:], label=f'Density at t={t5_hours} hours')
plt.plot(x, n_data[t6,:], label=f'Density at t={t6_hours} hours')
plt.plot(x, n_data[t7,:], label=f'Density at t={t7_hours} hours')
plt.plot(x, n_data[t8,:], label=f'Density at t={t8_hours} hours')
plt.plot(x, n_data[t9,:], label=f'Density at t={t9_hours} hours')

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Cell Density Spatial Profile')
plt.xlabel('x axis in m')
plt.ylabel(r'Cell Density in $\frac{\text{cell}}{\text{m}^3}$')
plt.tight_layout()
plt.show()


plt.figure()
x = np.linspace(-L/2, L/2, len(c_data[0,:]))

# plt.plot(x, c_data[0,:], label=f'EGF concen. at t={t0} hours')
# plt.plot(x, c_data[t1,:], label=f'EGF concen. at t={t1_hours} hours')
# plt.plot(x, c_data[t2,:], label=f'EGF concen. at t={t2_hours} hours')
# plt.plot(x, c_data[t3,:], label=f'EGF concen. at t={t3_hours} hours')
plt.plot(x, c_data[t4,:], label=f'EGF concen. at t={t4_hours} hours')
plt.plot(x, c_data[t5,:], label=f'EGF concen. at t={t5_hours} hours')
plt.plot(x, c_data[t6,:], label=f'EGF concen. at t={t6_hours} hours')
plt.plot(x, c_data[t7,:], label=f'EGF concen. at t={t7_hours} hours')
plt.plot(x, c_data[t8,:], label=f'EGF concen. at t={t8_hours} hours')
plt.plot(x, c_data[t9,:], label=f'EGF concen. at t={t9_hours} hours')

# Adding a title and labels
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('EGF Concen. Spatial Profile')
plt.xlabel('x axis in m')
plt.ylabel(r'EGF concentration in $\frac{\text{mol}}{\text{m}^3}$')
plt.tight_layout()
plt.show()



