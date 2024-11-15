# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 09:52:17 2024

@author: Kevin Roberts
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:39:54 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# EXAM REDO
# just redoing the problems that I missed on the exam


#############
# PROBELM 2 #
#############

# defining some parameters
C_e = 1 # (g/l)
C_f = 0.70*C_e # (g/l)
D = 0.33*10**(-4) # (m^2/day)
dx = 0.01
dy = 0.01
k_m = 0.93*10**(-3) 
Sh = k_m*dx/D
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo = D*dt/(dx**2)

save_step = 15
interval_steps = int(save_step / dt)

x = 0.2 # (cm)
y = 0.1 # (cm)

xn = int(x/dx)
yn = int(y/dy)

c = np.zeros((2,xn+1,yn+1)) # sets the initial conditions
data = np.zeros((1, xn+1, yn+1))
data[0] = c[1]

step_saved = 1

loop_counter = 0
time = 0
minimum = 0

while minimum < C_f:
    
    c[0,:,:] = c[1,:,:]   
    
    #interior
    for i in range(0,xn+1):
        for j in range(0,yn+1):
            if i == 0 and j == 0:  # upper left corner
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0 and j == yn:  # upper right corner
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == 0:  # bottom left corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == yn:  # bottom right corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0:  # top
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif j == 0:  # left side
                c[1,i,j] = 2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif i == xn:  # bottom
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            else:
                c[1,i,j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))

    
    loop_counter += 1
    time += dt
    
    minimum = np.min(c[1])
    
    # saving the data every 15 steps
    if loop_counter % interval_steps == 0:
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        minimum = np.min(c[1,:,:])
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")
        step_saved += 1

    # Check if all agar has at least target concentration
    if minimum >= C_f:
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        break

end_time = dt*loop_counter

print("QUESTION 2 ANSWER: " + str(end_time))
print()

#############
# PROBELM 4 #
#############

t_D = 0.2**2/(4*D)
print("QUESTION 4 ANSWER: " + str(t_D))
print()

#############
# PROBELM 6 #
#############

D = 0.33*(10**-4) # (m^2/day)
dx = 0.01 # (cm)
dy = 0.01 # (cm)
k_m = 0.93*10**(-3) 
C_e = 1 # (g/l)
C_f = 0.65*C_e # (g/l)
Sh = round(k_m*dx/D,3)
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo = D*dt/(dx**2)

save_step = 15
interval_steps = int(save_step / dt)

# defining the concentration
final_x = 0.2 # (m)
final_y = 0.1 # (m)

xn = int(final_x/dx)
yn = int(final_y/dy)

c = np.zeros((2,xn+1,yn+1)) # sets the initial conditions
data = np.zeros((1, xn+1, yn+1))
data[0] = c[0]

step_saved = 0

# adding the glucose on the side 

loop_counter = 0
time = 0
end_time = 0
minimum = 0

while minimum < C_f:
    
    c[0,:,:] = c[1,:,:]   
    
    #interior
    for i in range(0,xn+1):
        for j in range(0,yn+1):
            if i == 0 and j == 0:  # upper left corner
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0 and j == yn:  # upper right corner
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == 0:  # bottom left corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == yn:  # bottom right corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0:  # top
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif j == 0:  # left side
                c[1,i,j] = 2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif i == xn:  # bottom
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            else:
                c[1,i,j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))

    loop_counter += 1
    time += dt
    
    # saving the data every 15 steps
    if loop_counter % interval_steps == 0:
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        minimum = np.min(c[1,:,:])
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")
        

    # Check if all agar has at least target concentration
    if np.all(c[1,:,:] >= C_f):
        # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        end_time = dt*loop_counter
        break

end_time = dt*loop_counter

print("QUESTION 5 ANSWER: " + str(end_time))
print()
print("QUESTION 6 ANSWER: " + str(c[1,int(xn/2),int(yn/2)]*100))
print()

#############
# PROBELM 7 #
#############

D = 0.33*(10**-4) # (m^2/day)
dx = 0.01 # (cm)
dy = 0.01 # (cm)
k_m = 0.93*10**(-3) 
C_e = 1 # (g/l)
C_f = 0.65*C_e # (g/l)
Sh = round(k_m*dx/D,3)
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo = D*dt/(dx**2)

save_step = 15
interval_steps = int(save_step / dt)

# defining the concentration
final_x = 0.2 # (m)
final_y = 0.1 # (m)

xn = int(final_x/dx)
yn = int(final_y/dy)

c = np.zeros((2,xn+1,yn+1)) # sets the initial conditions
data = np.zeros((1, xn+1, yn+1))
data[0] = c[0]

step_saved = 0

# adding the glucose on the side 

loop_counter = 0
time = 0
end_time = 0
minimum = 0

while minimum < C_f:
    
    c[0,:,:] = c[1,:,:]   
    
    #interior
    for i in range(0,xn+1):
        for j in range(0,yn+1):
            if i == 0 and j == 0:  # upper left corner
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 1*2*Sh*C_e) + (1 - 4*Fo - 1*4*Fo*Sh)*c[0,i,j]
            elif i == 0 and j == yn:  # upper right corner
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == 0:  # bottom left corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 1*2*Sh*C_e) + (1 - 4*Fo - 1*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == yn:  # bottom right corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0:  # top
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == 0:  # left side
                c[1,i,j] = 2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif i == xn:  # bottom
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            else:
                c[1,i,j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))

    loop_counter += 1
    time += dt
    
    # saving the data every 15 steps
    if loop_counter % interval_steps == 0:
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        minimum = np.min(c[1,:,:])
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")
        

    # Check if all agar has at least target concentration
    if np.all(c[1,:,:] >= C_f):
        # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        end_time = dt*loop_counter
        break

end_time = dt*loop_counter

print("QUESTION 7 ANSWER: " + str(end_time))
print()

#############
# PROBELM 8 #
#############

print("QUESTION 8 ANSWER: " + str(c[1,int(xn/2),int(yn/2)]*100))
print()

##############
# PROBELM 11 #
##############

V = 1
A = 6*0.1**2
j = 0.3

C_final = (11 * V - j*A*60)/V
print("QUESTION 11 ANSWER: " + str(C_final))
print()


##############
# PROBELM 12 #
##############

D = 0.33*(10**-4) # (m^2/day)
dx = 0.01 # (cm)
dy = 0.01 # (cm)
k_m = 0.93*10**(-3) 
C_e = 1 # (g/l)
C_f = 0.65*C_e # (g/l)
Sh = round(k_m*dx/D,3)
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo = D*dt/(dx**2)

save_step = 15
interval_steps = int(save_step / dt)

# defining the concentration
final_x = 0.2 # (m)
final_y = 0.1 # (m)

xn = int(final_x/dx)
yn = int(final_y/dy)

c = np.zeros((2,xn+1,yn+1)) # sets the initial conditions
data = np.zeros((1, xn+1, yn+1))
data[0] = c[0]

step_saved = 0

# adding the glucose on the side 

loop_counter = 0
time = 0
end_time = 0
minimum = 0

while minimum < C_f:
    
    c[0,:,:] = c[1,:,:]   
    
    #interior
    for i in range(0,xn+1):
        for j in range(0,yn+1):
            if i == 0 and j == 0:  # upper left corner
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 1*2*Sh*C_e) + (1 - 4*Fo - 1*4*Fo*Sh)*c[0,i,j]
            elif i == 0 and j == yn:  # upper right corner
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == 0:  # bottom left corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 1*2*Sh*C_e) + (1 - 4*Fo - 1*4*Fo*Sh)*c[0,i,j]
            elif i == xn and j == yn:  # bottom right corner
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 0.5*2*Sh*C_e) + (1 - 4*Fo - 0.5*4*Fo*Sh)*c[0,i,j]
            elif i == 0:  # top
                c[1,i,j] = 2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == 0:  # left side
                c[1,i,j] = 2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif i == xn:  # bottom
                c[1,i,j] = 2*Fo*(c[0,i-1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j]
            else:
                c[1,i,j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))

    loop_counter += 1
    time += dt
    
    # saving the data every 15 steps
    if loop_counter % interval_steps == 0:
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        minimum = np.min(c[1,:,:])
        # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")
        

    # Check if all agar has at least target concentration
    if np.all(c[1,:,:] >= C_f):
        # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        end_time = dt*loop_counter
        break

end_time = dt*loop_counter
