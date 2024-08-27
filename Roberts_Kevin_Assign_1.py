"""
Created on Mon Aug 26 15:02:27 2024

@author: Kevin Roberts
"""

# answers to the homework problem questions:
    
# a. Where is mis-specification of dt most likely to cause error?
# 
# We discussed in class that mis-specifcation may cause error at the middle or 
# the end of all time t, however, it is most likely to cause error in the 
# beginning, when dt, because these are where the changes are happening most rapidly. 
#
# b. WHere does [RL](t) change most rapidly?
#
# In the beginning of the simulation, when time is small. The rapid changes 
# occur in the first few time steps.
#
# c. 
#
#
#
#
#
#


import math
import numpy as np
import matplotlib.pyplot as plt


# approximating the differential equation

kr = 1/500 # 1/seconds
delt = 0.1 # time interval, seconds
n = int(100/delt)+1 # total intervals
array = np.zeros((n, 3)) # define array
array[0, 2] = 100 # initialize array
time = 0
for i in range(1, n):
    array[i, 0] = i
    time += delt
    array[i,1] += time
    array[i, 2] = array[i-1,2] * (1 - delt*kr)


# setting array of values for the exact solution

i_c = 100 # this is the initial condition

def f(t):
    return i_c * math.exp(-kr*t)


exact_array = np.zeros(n)
for i in range(0, n):
    eval_at = i*delt
    exact_array[i] = f(eval_at) 

print(exact_array)
print(array)

x = np.arange(0, 100.1, 0.1)

approx_array = [row[2] for row in array] # gets the third column of the array

plt.plot(x, exact_array, label='Exact Solution')
plt.plot(x, approx_array, label='Approximated Solution')

# Adding a title and labels
plt.title('Approximated Solution vs Exact Solution')
plt.xlabel('Time in Seconds')
plt.ylabel('RL')

# Adding a legend to distinguish between the arrays
plt.legend()

# Display the plot
plt.show()





