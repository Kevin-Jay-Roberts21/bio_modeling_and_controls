# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 09:52:17 2024

@author: Kevin Roberts
"""

import numpy as np

# EXAM REDO
# just redoing the problems that I missed on the exam

#############
# PROBLEM 2 #
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

xn = int(x/dx) + 1
yn = int(y/dy) + 1

c = np.zeros((2,xn,yn)) # sets the initial conditions
data = np.zeros((1, xn, yn))
data[0] = c[1]

step_saved = 1

loop_counter = 0
time = 0
minimum = 0

while minimum < C_f:
    
    c[0,:,:] = c[1,:,:]   
    
    # interior
    for i in range(1, xn-1):
        for j in range(1, yn-1):
            c[1,i,j] = c[0,i,j] + Fo*(c[0,i+1,j] - 2*c[0,i,j] + c[0,i-1,j] + c[0,i,j+1] - 2*c[0,i,j] + c[0,i,j-1])

    # left 
    for j in range(1, yn-1):
        c[1,0,j] = 2*Fo*(c[0,1,j] + 0.5*(c[0,0,j-1] + c[0,0,j+1]) + Sh*C_e) + (1-4*Fo - 2*Fo*Sh)*c[0,0,j]

    # Right side
    for j in range(1, yn-1):
        c[1,xn-1,j] = 2*Fo*(c[0,xn-2,j] + 0.5*(c[0,xn-1,j-1] + c[0,xn-1,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,xn-1,j]

    # bottom
    for i in range(1, xn-1):
        c[1,i,0] = 2*Fo*(c[0,i,1] + 0.5*(c[0,i-1,0] + c[0,i+1,0]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,0]

    # top
    for i in range(1, xn-1):
        c[1,i,yn-1] = 2*Fo*(c[0,i,yn-2] + 0.5*(c[0,i-1,yn-1] + c[0,i+1,yn-1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,yn-1]

    # corners
    c[1,0,0] = 2*Fo*(c[0,1,0] + c[0,0,1] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,0,0]
    c[1,0,yn-1] = 2*Fo*(c[0,1,yn-1] + c[0,0,yn-2] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,0,yn-1]
    c[1,xn-1,0] = 2*Fo*(c[0,xn-2,0] + c[0,xn-1,1] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,xn-1,0]
    c[1,xn-1,yn-1] = 2*Fo*(c[0,xn-2,yn-1] + c[0,xn-1,yn-2] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,xn-1,yn-1]

    loop_counter += 1

    # Check if the minimum concentration in the slab reaches the target concentration
    if np.min(c[1, :, :]) >= C_f:
        break

end_time = dt*loop_counter

print("QUESTION 2 ANSWER: " + str(end_time))
print()

#############
# PROBLEM 4 #
#############

t_D = 0.1**2/(4*D)
print("QUESTION 4 ANSWER: " + str(t_D))
print()

#############
# PROBLEM 5 #
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
        # print(f"Day {int{time}}: Minimum m concentration = {minimum:.4f} g/L")
        

    # Check if all c has at least target concentration
    if np.all(c[1,:,:] >= C_f):
        # print(f"All m got to at least {C_f} g/L at day {time:.2f}")
        step_saved += 1
        data = np.append(data, np.zeros((1, xn+1, yn+1)), axis=0)
        data[step_saved] = c[1,:,:]
        end_time = dt*loop_counter
        break

end_time = dt*loop_counter

print("QUESTION 5 ANSWER: " + str(end_time))
print()

#############
# PROBLEM 6 #
#############

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
C_f = 0.35*C_e # (g/l)
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
                c[1,i,j] = C_e 
            elif j == 0:  # left side
                c[1,i,j] = C_e 
            elif j == yn:  # right side
                c[1,i,j] = 2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + 0*Sh*C_e) + (1 - 4*Fo - 0*2*Fo*Sh)*c[0,i,j]
            elif i == xn:  # bottom
                c[1,i,j] = C_e
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
        # print(f"Day {int{time}}: Minimum m concentration = {minimum:.4f} g/L")
        

    # Check if all c has at least target concentration
    if np.all(c[1,:,:] >= C_f):
        break

end_time = dt*(loop_counter-1)

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

C_e = 20  
C_f = 0.7*C_e
D = 0.33 * 10 ** -4  
k_m = 0.93 * 10 ** -3  

dx = 0.01
dy = 0.01  

# Slab dimensions
x = 0.2  
y = 0.1
z = 0.1  

# Grid size
xn = int(x/dx)
yn = int(y/dy)

# Sherwood number and Biot number
Sh = k_m * dx / D
Bi = Sh

# Stability criteria
dt_interior = (dx ** 2) / (4 * D)
dt_side = (dx ** 2) / (2 * D * (2 + Bi))
dt_corner = (dx ** 2) / (4 * D * (1 + Bi))
dt_max = min(dt_interior, dt_side, dt_corner)
dt = round(dt_max, 3)
Fo = D * dt / (dx ** 2)

# Initialize concentration array
c = np.zeros([2, xn + 1, yn + 1])
m = np.zeros([2, xn + 1, yn + 1])
   
# Simulation variables
loop_counter = 0  # Time steps
m_slab = 0

# m variables
V = 5*(1/1000)  # Convert to liters
C_i = V*C_e


# Volume of a single cell in the grid
v = dx*dy*z  # Assumes each cell has a depth of Lz

# Simulation loop
while np.min(c[1, :, :]) < C_f:

    c[0, :, :] = c[1, :, :]
    C_updated = (C_i - m_slab) / V
    C_e = C_updated
    
    for i in range(0, xn+1):
        for j in range(0, yn+1):
            if i == 0 and j == 0:  # upper left 
                c[1,i,j] = (2*Fo*(c[0,i+1,j] + c[0,i,j+1] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/4
            elif i == 0 and j == yn:  # upper right 
                c[1,i,j] = (2*Fo*(c[0,i,j-1] + c[0,i+1,j] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/4
            elif i == xn and j == 0:  # bottom left 
                c[1,i,j] = (2*Fo*(c[0,i-1,j] + c[0,i,j+1] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/4
            elif i == xn and j == yn:  # bottom right 
                c[1,i,j] = (2*Fo*(c[0,i-1,j] + c[0,i,j-1] + 2*Sh*C_e) + (1 - 4*Fo - 4*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/4
            elif i == 0:  # top
                c[1,i,j] = (2*Fo*(c[0,i+1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/2
            elif j == 0:  # left 
                c[1,i,j] = (2*Fo*(c[0,i,j+1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/2
            elif j == yn:  # right 
                c[1,i,j] = (2*Fo*(c[0,i,j-1] + 0.5*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/2
            elif i == xn:  # bottom
                c[1,i,j] = (2*Fo*(c[0,i-1,j] + 0.5*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo - 2*Fo*Sh)*c[0,i,j])
                m[1,i,j] = c[1,i,j]*v/2
            else:
                c[1,i,j] = (c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2)))
                m[1,i,j] = c[1,i,j]*v

    m_slab = np.sum(m[1, :, :])  

    loop_counter += 1

    if np.min(c[1, :, :]) >= C_f:
        break

end_time = loop_counter*dt

print("QUESTION 12 ANSWER: " + str(end_time))