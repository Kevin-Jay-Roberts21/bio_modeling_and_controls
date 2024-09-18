#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:07:38 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt

#############
# Problem 1 #
#############

# a) Determining dt
print("PROBLEM 1a OUTPUT:")

dt = 0.26 # picking an arbitrary dt
x_0 = 10 # initial prey
y_0 = 10 # initial predator
a = 1.1 # maximum growth rate per capita
b = 0.4 # effect of predator on prey death rate
g = 0.4 # predator death rate per capita
d = 0.1 # effect of prey on predator growth rate

j = 2
pred_perc_diff = 25 # initialization of perc diff, just some arbitrary value greater than 0.25
prey_perc_diff = 25 # initialization of perc diff, just some arbitrary value greater than 0.25

# starting a while loop to determine dt (using 25% as the percentage difference condition)
while (prey_perc_diff >= 25) or (pred_perc_diff >= 25):
    
    # calculating the next population value
    new_x = x_0 + dt*(a*x_0 - b*x_0*y_0)
    new_y = y_0 + dt*(-g*y_0 + d*x_0*y_0)
    
    # calculating the perc difference of x_i and x_(i-1)
    prey_perc_diff = 2*np.abs(new_x - x_0)/(new_x + x_0)*100
    
    # calculating the perc difference of y_i and y_(i-1)
    pred_perc_diff = 2*np.abs(new_y - y_0)/(new_y + y_0)*100
    
    # printing to the console dt for the first 5 loops
    if (j >= 2) and (j <= 7):
        print("The value of dt at iteration " + str(j-1) + " is " + str(dt))
        print("Perc diff of x and y at iteration " + str(j-1) + " is " + str(np.round(prey_perc_diff, 2)) + "% and " + str(np.round(pred_perc_diff, 2)) + "%, respectively.")

    # decreasing dt by 0.01 for each loop
    dt = np.round(dt - 0.01, 2)
    j += 1
 
final_dt = 0.06
      
print("The final dt value is: " + str(final_dt))
print()
print()

# b) Solving for x(t) and y(t) discretely

final_t = 100 # (years)

def get_pred_prey_data(time, dt):
    n = int(time/dt) + 1 # iternation number
    pred_prey_array = np.zeros((n, 4)) # define array (3rd and 4th col are prey and pred respectively)
    pred_prey_array[0, 2] = x_0
    pred_prey_array[0, 3] = y_0
    time = 0
    # solving for x and y discretely
    for i in range(1, n):
        pred_prey_array[i, 0] = i
        time += dt
        pred_prey_array[i, 1] += time
        pred_prey_array[i, 2] = pred_prey_array[i-1, 2] + dt*(a*pred_prey_array[i-1, 2] - b*pred_prey_array[i-1, 2]*pred_prey_array[i-1, 3])
        pred_prey_array[i, 3] = pred_prey_array[i-1, 3] + dt*(-g*pred_prey_array[i-1, 3] + d*pred_prey_array[i-1, 2]*pred_prey_array[i-1, 3])
        
    return pred_prey_array

array_dt_0 = get_pred_prey_data(final_t, final_dt)

# plotting
x = np.arange(0, len(array_dt_0)*final_dt, final_dt)

plt.plot(x, [row[2] for row in array_dt_0], label='Prey Pop.')
plt.plot(x, [row[3] for row in array_dt_0], label='Predator Pop.')

# Adding a title and labels
plt.legend()
plt.title('Predator and Prey Dynamics')
plt.xlabel('Time in Years')
plt.ylabel('Population')


# c) Examine the sensitivity of y(t) and x(t) to dt

# defining different dt values
dt_0 = np.round(final_dt*(3**0), 2) # already have this one
dt_1 = np.round(final_dt*(3**1), 2)
dt_2 = np.round(final_dt*(3**2), 2)
dt_3 = np.round(final_dt*(3**3), 2)
dt_4 = np.round(final_dt*(3**4), 2)

print(dt_0)
print(dt_1)
print(dt_2)
print(dt_3)
print(dt_4)






array_dt_1 = get_pred_prey_data(final_t, dt_1)
array_dt_2 = get_pred_prey_data(final_t, dt_2)
array_dt_3 = get_pred_prey_data(final_t, dt_3)
array_dt_4 = get_pred_prey_data(final_t, dt_4)


print(array_dt_2)



# plotting 4 subplots
x_dt_1 = np.arange(0, len(array_dt_1)*dt_1, dt_1)
x_dt_2 = np.arange(0, len(array_dt_2)*dt_2, dt_2)
x_dt_3 = np.arange(0, len(array_dt_3)*dt_3, dt_3)
x_dt_4 = np.arange(0, len(array_dt_4)*dt_4, dt_4)

plt.figure()

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
fig.suptitle("Plotting Predator and Prey Populations with Different dt")
ax1.plot(x_dt_1, [row[2] for row in array_dt_1], label="Prey Population") #
ax1.plot(x_dt_1, [row[3] for row in array_dt_1], label="Predator Population")
ax2.plot(x_dt_2, [row[2] for row in array_dt_2], label="Prey Population") #
ax2.plot(x_dt_2, [row[3] for row in array_dt_2], label="Predator Population")
ax3.plot(x_dt_3, [row[2] for row in array_dt_3], label="Prey Population") #
ax3.plot(x_dt_3, [row[3] for row in array_dt_3], label="Predator Population")
ax4.plot(x_dt_4, [row[2] for row in array_dt_4], label="Prey Population") #
ax4.plot(x_dt_4, [row[3] for row in array_dt_4], label="Predator Population")
ax1.set_xlabel('Time in Years')
ax1.set_ylabel('Poplation Count')
ax2.set_xlabel('Time in Years')
ax2.set_ylabel('Poplation Count')
ax3.set_xlabel('Time in Years')
ax3.set_ylabel('Poplation Count')
ax4.set_xlabel('Time in Years')
ax4.set_ylabel('Poplation Count')
ax1.set_title("dt = " + str(np.round(dt_1, 2)))
ax2.set_title("dt = " + str(np.round(dt_2, 2)))
ax3.set_title("dt = " + str(np.round(dt_3, 2)))
ax4.set_title("dt = " + str(np.round(dt_4, 2)))
plt.subplots_adjust(hspace=0.9)
plt.legend()
plt.tight_layout()
plt.show()

# NOTE: The graphs may not seem to make sense in langer time shots. To
# see the interactions for the different dt's, it may be useful to decrease
# the total time to a smaller value.



# d) Evaluate the sensitivity of y(t) and x(t) to dt

x = np.array([dt_0, dt_1, dt_2, dt_3, dt_4])

dt_0_x_max = max([row[2] for row in array_dt_0])
dt_1_x_max = max([row[2] for row in array_dt_1])
dt_2_x_max = max([row[2] for row in array_dt_2])
dt_3_x_max = max([row[2] for row in array_dt_3])
dt_4_x_max = max([row[2] for row in array_dt_4])
dt_0_y_max = max([row[3] for row in array_dt_0])
dt_1_y_max = max([row[3] for row in array_dt_1])
dt_2_y_max = max([row[3] for row in array_dt_2])
dt_3_y_max = max([row[3] for row in array_dt_3])
dt_4_y_max = max([row[3] for row in array_dt_4])






largest_x_array = np.array([dt_0_x_max, dt_1_x_max, dt_2_x_max, dt_3_x_max, dt_4_x_max])
largest_y_array = np.array([dt_0_y_max, dt_1_y_max, dt_2_y_max, dt_3_y_max, dt_4_y_max])
                    
print(largest_x_array)
print(largest_y_array)


                                                                                                                          
fig, ax1 = plt.subplots()
fig.suptitle("Max of x and y for each dt.")
line1 = ax1.scatter(x, largest_x_array, c="blue", marker="D", label="max(x(t))") # Plot (x^3 vs x)
ax2 = ax1.twinx()
line2 = ax2.scatter(x, largest_y_array, c="red", linestyle="solid", label="max(y(t))")
ax1.set_xlabel('dt value')
ax1.set_ylabel('highest x(t) value', color = 'blue')
ax2.set_ylabel('highest y(t) value', color = 'red')
ax1.tick_params('y', colors = 'blue')
ax2.tick_params('y', colors = 'red')
ax1.legend((line1,line2),('max(x(t))', 'max(y(t))'), loc = 'upper right')
ax1.set_ylim(-100,5000)
ax2.set_ylim(-100,5000)
plt.tight_layout()
plt.show()






#############
# Problem 2 #
#############

# a) Determine dt

# first I'll deine all the parameters being used
alpha_max = 50*10**(-6) # (M/s)
alpha_min = 1*10**(-6) # (M/s)
Ca = 0.66*10**(-6) # (M)
K_cyc = 0.135*10**(-6) # (M)
m_cyc = 2 # (dimensionless)
alpha = alpha_min + (alpha_max-alpha_min)/(1 + (Ca/K_cyc)**(m_cyc)) # (M/s)
V_cyt = 1076 # (um^3)
A_disc = V_cyt/0.007 # (um)
k_hyd = 7*10**(-5) # (um^3/s)
k_hyd_plus = 1*10**(-5) # (um^3/s)
PDE = 100 # (um^(-2))
J_dark = 66 #(pA)
A_activ = 800*np.pi*(5.5**(2)) # (um^2)
v_RE = 195 # (1/s)
k_R = 2.6 # (1/s)
k_E = 0.6 # (1/s)
phi = 1 # (dimensionless)
J_cGMP_max = 7*10**(3) # (pA)
k_cGMP = 32*10**(-6) # (M)
mcGMP = 2 # (dimensionless)

cGMP_0 = 3*10**(-6) # (M)

def E_plus(t):
    E_at_t = phi*v_RE/(k_R - k_E)*(np.e**(-k_E*t) - np.e**(-k_R*t))        
    return E_at_t

def PDE_plus(t):
    PDE_plus_at_t = E_plus(t)/(2*A_activ) 
    return PDE_plus_at_t


j = 2
initial_dt = 0.2

# starting a while loop to determine dt
while True:
    
    # calculating the next population value
    new_cGMP = cGMP_0 + initial_dt*(alpha - A_disc/V_cyt*(k_hyd*PDE - k_hyd_plus*PDE_plus(0))*cGMP_0)
        
    if new_cGMP >= cGMP_0*0.999:
        cGMP_final_dt = initial_dt
        break
    
    # decreasing dt by 0.01 for each loop
    initial_dt = np.round(initial_dt - 0.01, 2)
    j += 1
    
    # printing to the console dt for the first 5 loops
    if (j >= 2) and (j <= 7):
        print("The value of dt at iteration " + str(j-1) + " is " + str(initial_dt))
        print("The new value of cGMP at " + str(j-1) + " is " + str(new_cGMP))
 
      
print("The final dt value is: " + str(cGMP_final_dt))
print()
print()


# b) Solve for [cGMP](t) using equation 3










# c) Examine parametric sensitivity of [cGMP](t) to k_R and k_E in terms of precision and response time
