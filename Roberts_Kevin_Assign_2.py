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
k_f = (10**5) # (1/(M*s)) the forward rate constant
r = 5*(10**4) # (1/(M*s)) the ratio k_cat/k_m
k_m = k_cat/r # (M) since we know that k_cat/k_m = 5*10^4 (1/(M*s))
k_r = k_cat*k_f/r - k_cat # (1/s) using the fact that k_m = (k_r+k_cat)/k_f

E_out = 70*(10**(-9)) # (M) molars
S_initial = 200*(10**(-9)) # (M) molars
V_i = 70 # (microliters)

# verifying the (i) condition
def check_i_condition(E_o):
    
    check = False
    
    # defining the LHS and RHS of (i)
    i_LHS = E_o/(k_m + S_initial) 
    i_RHS = (1 + k_r/k_cat)*(1 + S_initial/k_m)
    
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
    ii_LHS = E_o/(k_m + S_intial)
    ii_RHS = 1    
    
    if ii_LHS < ii_RHS*0.01:
        check = True
        return check
    else:
        return check

def get_E_o_and_volume(E_o):
    if check_i_condition(E_o) and check_i_condition(E_o):
        print("Both (i) and (ii) are satisfied.")
        
        # we need to solve the volume V_o of the eqn: E_o*(V_i + V_o) = E_out*V_o
        # which is: V_o = -E_o*V_i/(E_o - E_out)
        volume = -E_o*V_i/(E_o - E_out)
        
        print("The value of E_o is: " + str(E_o) + " M and the volume is: " + str(round(volume, 2)) + " microliters.")
    else:
        print("Either (i) or (ii) or both was not satisfied.")
        
    return volume    
        
enzyme_to_test = 22*(10**(-9)) # (M)
E_0 = enzyme_to_test
V_add = get_E_o_and_volume(enzyme_to_test)



# Finding S_0 (iterating through to find S_0)
S_1 = S_initial*V_i/(V_add + V_i)
S_2 = S_1*V_i/(V_add + V_i)
S_3 = S_2*V_i/(V_add + V_i)
S_4 = S_3*V_i/(V_add + V_i)
S_5 = S_4*V_i/(V_add + V_i)
S_6 = S_5*V_i/(V_add + V_i)
print(S_1)
print(S_2)
print(S_3)
print(S_4)
print(S_5)
print(S_6)
S_0 = S_6

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
# d. If the quassi steady state of ES is ES = E_0*S_0/(k_m + S_0), then when 
#    S_0 changes little (say we add a few nM of concentration), then we end up seeing 
#    that ES tends to approach a steady value. I've printed out my results under 
#    PROBLEM 2 OUTPUT
# 
#
# e. Just like the last homework assignment, a reasonable increment in time is 
#    one that does not change the previous value from the new value very much. We 
#    typically want the change to vary from the old to the new value about 10%. So 
#    in this situation we choose the dt to be dt = 0.01 seconds
# 
# 
# f. To find how many seconds in this problem, we must solve for t in the following:
#    [ES]/[ES]_qausi = 1 - exp(-t/t_qssa). Solving for t results in the following:
#    t = -t_qssa*ln(0.01) = , when [ES]/[ES]_qausi = 0.01. I printed out the results.
# 
# 
# g. The value of n that will achieve 99% of quasi-steady state is the following:
#    n = (final_t/dt) + 1 = 

ES_quasi = E_0*S_0/(k_m + S_0)
t_qssa = k_f*(k_m + S_0)
dt = 0.01
print("The following is value of ES steady state when S increases by a nM, 2nM, and 3nM:")
print("Original ES: " + str(ES_quasi))
print("ES plus one nM: " + str(E_0*(S_0 + (10**(-9)))/(k_m + S_0 + (10**(-9)))))
print("ES plus two nM: " + str(E_0*(S_0 + (20**(-9)))/(k_m + S_0 + (20**(-9)))))
print("ES plus three nM: " + str(E_0*(S_0 + (30**(-9)))/(k_m + S_0 + (30**(-9)))))
print()
print()
print("The number of seconds it takes for the transition to reach quasi-steady-state 99% completness:")
final_t = -t_qssa*np.log(0.01)
print("t = " + str(final_t))
print()
n = int(final_t/dt) + 1
print("The value of n that will achieve 99% of the quasi-steady state: n = " + str(n))
print()
print()








#############
# Problem 3 #
#############
print("PROBLEM 3 OUTPUT:")

# defining delta_t and n

# approximating the differential equation with dt = 0.1
E_0 = enzyme_to_test # (M) this is the value satisfied by criteria (ii)
ES_0 = 0 # (M) (this is an arbitrary value)
ES_array = np.zeros((n, 4)) # define array (4th row will be [ES]_qssa)
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
plt.title('Comparing ES fin. diff with ES Mich. Men.')
plt.xlabel('Time in Seconds')
plt.ylabel('ES in Molars')

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




print()
print()

#############
# Problem 6 #
#############
print("PROBLEM 6 OUTPUT:")
print("Told to skip question 6")






