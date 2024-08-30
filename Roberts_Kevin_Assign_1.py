"""
Created on Mon Aug 26 15:02:27 2024

@author: Kevin Roberts
"""

# answers to the homework problem questions:
    
# a.
# We discussed in class that mis-specifcation may cause error at the middle or 
# the end of all time t, however, it is most likely to cause error in the 
# beginning, when dt, because these are where the changes are happening most rapidly. 
#
# b.
# In the beginning of the simulation, when time is small. The rapid changes 
# occur in the first few time steps. This is for the exact solution, and for the
# other two solutions.
#
# c. 
# Below, I've displayed the rate of change for the exact solution, the approximated
# solution with dt = 0.1s and the approximated solution for dt = 20s. For the 
# function y = 100exp(-3/500t), the greatest change happens at time 0. Taking the
# derivative of y, and plugging in 0, we find the rate of change to be 0.60
#
# d.
# The rate of change here comes from the second derivative of the y function 
# described below, I found this to be: 0.0036
#
# e. A reasonable increment in time is one that does not change the previous value
# from the new value very much. As discussed in class, we typically want the change
# to vary from the old to the new value about 10%.
# 
# f.
# The value of n that will reduce the receptor - ligand complex concentration 
# by 91% results from solving t in the following equation: 9.1 = 100*exp(-3/500 * t)
# which is approximated: f_t = -500/3 * ln(9.1/100) = 399.48 secs. But we are looking
# for n = (f_t/dt)+1 = (399.48/0.1)+1 = 3995.8 

import math
import numpy as np
import matplotlib.pyplot as plt



kr = 3/500 # 1/seconds
i_c = 100 # defining the initial condition

# approximating the differential equation with dt = 0.1
delt_1 = 0.1 # time interval, seconds
n_1 = int(100/delt_1)+1 # total intervals
array_1 = np.zeros((n_1, 3)) # define array
array_1[0, 2] = i_c # initialize array
time = 0
for i in range(1, n_1):
    array_1[i, 0] = i
    time += delt_1
    array_1[i,1] += time
    array_1[i, 2] = array_1[i-1,2] * (1 - delt_1*kr)

# approximating the differential equation with dt = 2
delt_2 = 20 # time interval, seconds
n_2 = int(100/delt_2)+1 # total intervals
array_2 = np.zeros((n_2, 3)) # define array
array_2[0, 2] = i_c # initialize array
time = 0
for i in range(1, n_2):
    array_2[i, 0] = i
    time += delt_2
    array_2[i,1] += time
    array_2[i, 2] = array_2[i-1,2] * (1 - delt_2*kr)


def f(t):
    return i_c * math.exp(-kr*t)

exact_array = np.zeros(n_1)
for i in range(0, n_1):
    eval_at = i*delt_1
    exact_array[i] = f(eval_at) 

x_1 = np.arange(0, n_1*delt_1, delt_1)
x_2 = np.arange(0, n_2*delt_2, delt_2)
approx_array_1 = [row[2] for row in array_1] # gets the third column of the array
approx_array_2 = [row[2] for row in array_2] # gets the third column of the array


plt.plot(x_1, exact_array, label='Exact Solution')
plt.plot(x_1, approx_array_1, label='Approx Solution dt=0.1')
plt.plot(x_2, approx_array_2, label='Approx Solution dt=20')

# Adding a title and labels
plt.title('Approximated Solutions vs Exact Solution')
plt.xlabel('Time in Seconds')
plt.ylabel('RL')

# Adding a legend to distinguish between the arrays
plt.legend()

# Display the plot
plt.show()


def get_most_rapid_change(array, step_size):
    
    greatest_change = 0
    location = [0, 0]
    for i in range(len(array)-1):
        if np.abs(array[i+1] - array[i]) > greatest_change:
            greatest_change = abs(array[i+1] - array[i])
            location[0] = i*step_size
            location[1] = i+1*step_size
            
    print("The most rapid change is: " + str(greatest_change) + " and it occurs from " + str(location[0]) + "s to " + str(location[1]) + "s." )

get_most_rapid_change(exact_array, delt_1)
get_most_rapid_change(approx_array_1, delt_1)
get_most_rapid_change(approx_array_2, delt_2)




