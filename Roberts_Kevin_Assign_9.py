#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:28:48 2024

@author: kevinjayroberts
"""
import numpy as np
import matplotlib.pyplot as plt

#############
# PROBLEM 1 #
#############

# a) Console Output

z = 0.4 # length (m)
y = 0.1 # thickness (m)
x = 0.1 # height (m)

D = 0.252*(10**-4) # (m^2/day)
dx = 0.01 # (cm)
dy = 0.01 # (cm)
k_m = 3.1536*10**(-3) 
C_e = 20 # (g/l)
C_f = 15 # (g/l)
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
final_x = 0.4 # (m)
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

# printing all of the output stuff
print("PROBELM 1a OUTPUT:")
print("x[m]: " + str(round(c[1][xn-1][0],3)))
print("y[m]: " + str(round(c[1][1][yn-1],3)))
print("Aspect Ratio: 40x10x10")
print("Diffusion time to >= 15 g/l: " + str(loop_counter*dt))
print("Sherwood Number: " + str(Sh))
for i in range(len(data)-1):
    minimum = np.min(data[i,:,:])
    time = i*15
    print(f"Day {int(time)}: Minimum analgesic concentration = {minimum:.3f} g/L")
end_min = np.min(data[len(data)-1,:,:])
print(f"Day {int(end_time)} (last day): Minimum analgesic concentration = {end_min:.3f} g/L")
print()
print()
# b) Figure 1a

# getting the all the minimum values at each 15 days
mins = []

for i in range(len(data)):
    mins.append(np.min(data[i]))

L = 0.1 # (m)
t = np.arange(0, (len(data)-1)*save_step, save_step)
t= np.append(t, end_time)
dim_less_time = [round(i * (D / L**2), 3) for i in t]

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# First subplot for the original time plot
ax1.plot(t, mins, marker='o', label='Analgesic minimums')

# Customizing the first subplot
ax1.set_title("Plotting Analgesic Mins at every 15 days")
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Minimum Analgesic values (g/l)")
ax1.set_xticks(t)  # Set x-ticks to match the time values
ax1.legend()
ax1.grid()

# Second subplot for the dimensionless time plot
ax2.plot(dim_less_time, mins, marker='o', label='Analgesic minimums')

# Customizing the second subplot
ax2.set_title("Plotting Analgesic Mins (Dimensionless Time)")
ax2.set_xlabel("Dimensionless Time")
ax2.set_ylabel("Minimum Analgesic values (g/l)")
ax2.tick_params(axis='x', rotation=45)
ax2.set_xticks(dim_less_time)  # Set x-ticks to match the time values
ax2.legend()
ax2.grid()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()

# c) Explain
print("PROBLEM 1c OUTPUT")
print("i) The concentration in assignment 8 tooka  lot longer to get to a minimum of 15 g/l than the concentration in assignment 9, even though the target concentration in assignment 8 was 20 g/l.")
print("ii) We could cahnge the diffusion constant, or we could also change the boundary conditions.")
print("iii) Changing the diffusion constant would make the concentration in assignment 9 diffuse slower to match assignment 8. And of course if we change the boundary conditions then assignment 8 would exactly match assignment 9.")
print()
print()

# d) Figure 1b 
# Plot concentration profiles in 3D
x_grid = np.linspace(0, final_x, xn+1)  # x-axis points
y_grid = np.linspace(0, final_y, yn+1)  # y-axis points
X, Y = np.meshgrid(x_grid, y_grid)      # create the meshgrid

for m in range(len(data)-1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface, with correctly matched X, Y, and data[m,:,:] shapes
    ax.plot_surface(X, Y, data[m,:,:].T, cmap='viridis')

    ax.set_title(f'Analgesic at Day {m * save_step}')
    ax.set_xlabel('Length in cm')
    ax.set_ylabel('Thickness in cm')
    ax.set_zlabel('Concentration in g/L')
    plt.show()
# plotting the last day
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface, with correctly matched X, Y, and data[m,:,:] shapes
ax.plot_surface(X, Y, data[len(data)-1,:,:].T, cmap='viridis')

ax.set_title(f'Analgesic at Day {int(end_time)}')
ax.set_xlabel('Length in cm')
ax.set_ylabel('Thickness in cm')
ax.set_zlabel('Concentration in g/L')
plt.show()

# e) Explain
print("PROBLEM 1e OUTPUT")
print("i) The concentration regardless of dimensional or dimensionless time or concentration, behaves similarly. It is important to note that the nice thing about nondimensionalizing a system is that we are easily able to model any situation and get similar results.")
print("ii) The axial concentration profiles in Figure 1b show a sharper decrease on one of the sides of the slab and this is due to the modified boundary condition on the 'right' side, and the 0.5 being multiplied to the sherwood number in the discretization. They are similar because they have the overall same diffusion discretization method being applied to the center of the slab, and have the same initial conditions.")
print("iii) We might say the implications since there is no net mass transfer along one of the axis, the diffusion along that axis will not need to be explicitly included in a model.")
