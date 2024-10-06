#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 14:24:18 2024

@author: kevinjayroberts
"""

import matplotlib.pyplot as plt
import numpy as np

#############
# PROBELM 1 #
#############

# a) Use equation (1) to model glucose diffusing into the agar.

# defining some parameters
length = 20 # (cm)
thickness = 10 # (cm)
height = 1 # (cm)
C = 30 # (g/l)
C_f = 20 # (g/l)
D = 0.252*10**(-4) # (m^2/day)
dx = 1 # (cm)
dy = 1 # (cm)
k_m = 3.1536*10**(-3) 
Sh = k_m*dx/D
interior = dx**2/(4*D)
side = dx**2/(2*(2 + Sh))
corner = dx**2/(4*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = max(dt_maxes)
dt = round(dt_max, 3)

Fo_m = D*dt/(dx**2)

print("PROBELM 1a OUTPUT")
print("The stability limit for the interior is " + str(interior))
print("The stability limit for the side is: " + str(side))
print("The stability limit for the corner node is: " + str(corner))
print("The value of dt is: " + str(dt))
print("The value of Fo_m is: " + str(Fo_m))
print("The value of Sh is: " + str(Sh))
print()
print()


# b) Use dt = round(dt, 3) from 1.a. Use an array with 2*x*y dimensions to calculate 
#    agar concentrations at each node at each successive time.

# defining the concentration
final_x = 20 # (cm)
final_y = 10 # (cm)

xn = int(final_x/dx)
yn = int(final_y/dy)

c = zeros((2,xn,yn))

loop_counter = 1
while max(c()):
    for i in range(xn):
        for j in range(yn):
            # corners
            if i == 0 and j == 0:
                pass
            elif i == 0 and j == yn-1:
                pass
            elif i == xn-1 and j == 0:
                pass
            elif i == xn-1 and j == yn-1:
                pass
            
            # left side
            elif i == 0:
                pass
            
            # right side
            elif i == xn-1:
                pass
            
            # bottom
            elif j == 0:
                c(1, i, j) = 2*Fo_m*(c(0,i+1,j) + 1/2*(c(0,i,j) + c(0,i,j+1)) + Sh*k_m) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c(0,i,j)
            
            # top
            elif j == yn-1:
                c(1, i, j) = 2*Fo_m*(c(0,i+1,j) + 1/2*(c(0,i,j) + c(0,i,j+1)) + Sh*k_m) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c(0,i,j)
            else:
                c(1, i, j) = c(0,i,j) + D*dt*((c(0,i-1,j) - 2*c(0,i,j) + c(0,i+1,j))/(dx**2) + (c(0,i,j-1) - 2*c(0,i,j) + c(0,i,j+1))/(dy**2))




# c) Console Output






