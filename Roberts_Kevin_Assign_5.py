#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:27:40 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt

#############
# PROBLEM 1 #
#############

# a) Simulate u(x,t) 

# initial conditions: u = 2x 0<=x<=1/2, u = 2(1-x) 1/2<=x<=1
# boundary conditions: u_x=0 = u_x=1 = 0

# defining some variables
dt = 1/1000
dx = 1/10
alpha = 1
L = 1

# space and time intervals
x_0 = 0
x_f = 1
t_0 = 0
t_f = 0.04










# b) Output the simulation

















#############
# PROBLEM 2 #
#############

# a) Calculate u(x,t)









# b) Output calculated u(x,t)





# c) Describe comparing simulation vs. exact solution















#############
# PROBLEM 3 #
#############

# a) Compare simulated u(x,t)















#############
# PROBLEM 4 #
#############

# a) Use equation 1 with initial and boundary conditions in Problem 1 to simulate
# u(x,t) for 0<=x<=1 at times for which u >= 0.55. Use dt_max. Inside a while loop,
# use a two-row array to compute successive values of u(x,t)










# b) Output simulated u(x,t)





# c) Console output





# d) Describe loop accuracy

