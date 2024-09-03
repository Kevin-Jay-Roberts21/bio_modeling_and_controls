#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:25:52 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt

#############
# Problem 1 #
#############
print("PROBLEM 1 OUTPUT:")

# NOTE: I'm multiplying some of the concentrations by 10^{-9} to convert everything to 
# molars

# defining the rates
k_cat = 0.1 # (1/s) the caralytic rate
k_f = (10**5) # (1/(nM*s)) the forward rate constant
a = 5*(10**4) # (1/(nM*s)) the ratio k_cat/k_m
k_m = k_cat/a # (M) since we know that k_cat/k_m = 5*10^4 (1/(M*s))
k_r = (k_cat*k_f - a*k_cat)/a # (1/s) using the fact that k_m = (k_r+k_cat)/k_f

E_out = 70*(10**(-9)) # (M) nanomolars
S_0 = 200*(10**(-9)) # (nM) nonmolars
V_s = 70 # (microliters)

# verifying the (i) condition
def check_i_condition(E_o):
    
    check = False
    
    # defining the LHS and RHS of (i)
    i_LHS = E_o/(k_m + S_0) 
    i_RHS = (1 + k_r/k_cat)*(1 + S_0/k_m)
    
    # multiplying the i_RHS by 0.01 because we are defining the LHS to be 
    # "much less than" if it is less than a hundreth of the RHS
    if i_LHS < i_RHS*0.01:
        check = True
        return check
    else:
        return check

# verifying the (ii) condition
def check_ii_condition(E_o):
    
    check = False
    
    # defining the LHS and RHS of (ii)
    ii_LHS = E_o/(k_m + S_0)
    ii_RHS = 1    
    
    if ii_LHS < ii_RHS*0.01:
        check = True
        return check
    else:
        return check

def get_E_o_and_volume(E_o):
    if check_i_condition(E_o) and check_i_condition(E_o):
        print("Both (i) and (ii) are satisfied.")
        
        # we need to solve the volume V_o of the eqn: E_o*(V_s + V_o) = E_out*V_o
        # which is: V_o = -E_o*V_s/(E_o - E_out)
        volume = -E_o*V_s/(E_o - E_out)
        
        print("The value of E_o is: " + str(E_o) + " nM and the volume is: " + str(round(volume, 2)) + " microliters.")
    else:
        print("Either (i) or (ii) or both was not satisfied.")
        
        
enzyme_to_test = 22*(10**(-9)) # (M)
get_E_o_and_volume(enzyme_to_test)
print()
print()



#############
# Problem 2 #
#############
print("PROBLEM 2 OUTPUT:")
# a. The mis-specification of dt is most likely to cause error near the beginning,
#    when time is close to 0, because that's where the most rapid changes occur.
# 
# b. [ES](t) changes most rapidly at the beginning, when time is close to 0.
# 
# c. To find the rate of change of [ES](t), we must take the derivative of the
#    equation and evaluate it at 0, which is: 
#    ans = k_f*(E_0 - ES(0))*S_0 - k_r*ES(0) - k_cat*ES(0) -->
#    ans = 0.0001*(70 - 0)*200 - 0.1*0 - 0.1*0 = 1.4 (assuming that ES is 0 when 
#    time t = 0) 
# 
# d. 
# 
#
# e. Just like the last homework assignment, a reasonable increment in time is 
# one that does not change the previous value from the new value very much. We 
# typically want the change to vary from the old to the new value about 10%. So 
# in this situation we choose the dt to be dt = 0.1 seconds
# 
# 
# f. 
# 
# 
# g. 
print()
print()








#############
# Problem 3 #
#############
print("PROBLEM 3 OUTPUT:")

# defining delta_t and n

# approximating the differential equation with dt = 0.1
dt = 0.1 # time interval, seconds
final_t = 9
n = int(final_t/dt)+1 # total intervals
E_0 = enzyme_to_test # (nM) this is the value satisfied by criteria (ii)
ES_0 = 0 # (nM) (this is an arbitrary value)
ES_array = np.zeros((n, 4)) # define array
ES_array[0, 2] = ES_0 # initialize array
time = 0

for i in range(1, n):
    ES_array[i, 0] = i
    time += dt
    ES_array[i, 1] += time
    ES_array[i, 2] = dt*(k_f*(E_0 - ES_array[i-1,2])*S_0 - k_r*ES_array[i-1,2] - k_cat*ES_array[i-1, 2]) + ES_array[i-1, 2]


ES_M_M_const = E_0*S_0/(k_m + S_0)

for i in range(0, n):
    ES_array[i, 3] = ES_M_M_const

x = np.arange(0, n*dt, dt)

plt.plot(x, [row[2] for row in ES_array], label='ES Finite Difference')
plt.plot(x, [row[3] for row in ES_array], label='ES Michaelis-Menten')

# Adding a title and labels
plt.title('Comparing ES finite difference with ES Approximation (3)')
plt.xlabel('Time in Seconds')
plt.ylabel('RL')

# adding some scaling
# lower_y = 0
# upper_y = np.max([row[2] for row in ES_array])
# lower_x = 0
# upper_x = final_t
# plt.ylim(lower_y, upper_y)
# plt.xlim(lower_x, upper_x)


# Adding a legend to distinguish between the arrays
plt.legend()

# Display the plot
plt.show()

print()
print()


#############
# Problem 4 #
#############
print("PROBLEM 4 OUTPUT:")

# printing out the data to the console
print("Data from the ES Finite Difference and ES Michaelis-Menten")
print("j, time, [ES]_Fin_Diff, [ES]_Mich_Ment")

np.set_printoptions(precision=3, threshold=100, edgeitems=5, suppress=False)

print(ES_array)
print()
print()



#############
# Problem 5 #
#############
print("PROBLEM 5 OUTPUT:")


t_qssa = k_f*(k_m + S_0)


print()
print()

#############
# Problem 6 #
#############
print("PROBLEM 6 OUTPUT:")






