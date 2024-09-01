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

# NOTE: I'm multiplying some of the rates by 10^{-9} to convert everything to 
# nanomolars

# defining the rates
k_cat = 0.1 # (1/s) the caralytic rate
k_f = (10**5)*(10**(-9)) # (1/(nM*s)) the forward rate constant
a = 5*(10**4)*(10**(-9)) # (1/(nM*s)) the ratio k_cat/k_m
k_m = k_cat/a # (nM) since we know that k_cat/k_m = 5*10^4 (1/(M*s))
k_r = (k_cat*k_f - a*k_cat)/a # (1/s) using the fact that k_m = (k_r+k_cat)/k_f

print(k_r)
print(k_m)
print(k_f)



E_out = 70 # (nM) nanomolars
S_0 = 200 # (nM) nonmolars
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
        
        
enzyme_to_test = 22 # (nM)
get_E_o_and_volume(enzyme_to_test)




#############
# Problem 2 #
#############

# a. The mis-specification of dt is most likely to cause error near the beginning,
#    when time is close to 0, because that's where the most rapid changes occur.
# 
# b. [ES](t) changes most rapidly at the beginning, when time is close to 0.
# 
# c. To find the rate of change of [ES](t), we must take the derivative of the
#    equation and evaluate it at 0, which is: 
#    ans = k_f*E(0)*S(0) - k_r*ES(0) - k_cat*ES(0) -->
#    ans = 0.0001*22*70 - 0.1*??? - 0.1*??? =  











# defining delta_t and n

# approximating the differential equation with dt = 0.1
dt = 0.1 # time interval, seconds
n = int(100/dt)+1 # total intervals
E_0 = enzyme_to_test # (nM)
ES_initial = 200 # (nM) (this is an arbitrary value)
P_initial = 10 # (nM) (this is an arbitrary value)
S_initial = S_0 # (nM) (this is an arbitrary value)

ES_array = np.zeros((n, 3)) # define array
P_array = np.zeros((n, 3))
E_array = np.zeros((n, 3))
S_array = np.zeros((n, 3))

ES_array[0, 2] = ES_initial # initialize array
P_array[0, 2] = P_initial
E_array[0, 2] = E_0
S_array[0, 2] = S_initial

time = 0

for i in range(1, n):
    ES_array[i, 0] = i
    P_array[i, 0] = i
    E_array[i, 0] = i    
    S_array[i, 0] = i
    
    time += dt
    ES_array[i, 1] += time
    P_array[i, 1] += time
    E_array[i, 1] += time    
    S_array[i, 1] += time
    
    ES_array[i, 2] = dt*(k_f*E_array[i-1, 2]*S_array[i-1, 2] - k_r*ES_array[i-1,2] - k_cat*ES_array[i-1, 2]) + ES_array[i-1, 2]
    P_array[i, 2] = dt*(k_cat*ES_array[i-1, 2]) + P_array[i-1, 2]
    E_array[i, 2] = E_0 - ES_array[i-1, 2]
    S_array[i, 2] = dt*(-k_f*E_array[i-1, 2]*S_array[i-1, 2] + k_r*ES_array[i-1, 2]) + S_array[i-1, 2]


ES_exact = np.zeros((n, 3))
def ES_exact_ftn(time_step):
    return E_0*S_array[i, 2]/(k_m + S_array[i, 2])


x = np.arange(0, n*dt, dt)

plt.plot(x, [row[2] for row in ES_array], label='ES')
plt.plot(x, [row[2] for row in P_array], label='P')
plt.plot(x, [row[2] for row in E_array], label='E')
plt.plot(x, [row[2] for row in S_array], label='S')

# Adding a title and labels
plt.title('Stuff')
plt.xlabel('Time in Seconds')
plt.ylabel('RL')

# Adding a legend to distinguish between the arrays
plt.legend()

# Display the plot
plt.show()




