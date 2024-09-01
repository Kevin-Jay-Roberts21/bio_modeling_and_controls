#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:25:52 2024

@author: kevinjayroberts
"""
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

E_0 = 70 # (nM) nanomolars
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
        
        # we need to solve the volume V_o of the eqn: E_o*(V_s + V_o) = E_0*V_o
        # which is: V_o = -E_o*V_s/(E_o - E_0)
        volume = -E_o*V_s/(E_o - E_0)
        
        print("The value of E_o is: " + str(E_o) + " nM and the volume is: " + str(round(volume, 2)) + " microliters.")
    else:
        print("Either (i) or (ii) or both was not satisfied.")
        
        
enzyme_to_test = 22 # (nM)
get_E_o_and_volume(enzyme_to_test)




#############
# Problem 2 #
#############




