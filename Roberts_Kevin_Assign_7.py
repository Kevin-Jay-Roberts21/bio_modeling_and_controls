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
C_e = 30 # (g/l)
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

dx =1
dy = 1


# defining the concentration
final_x = 20 # (cm)
final_y = 10 # (cm)

xn = int(final_x/dx)
yn = int(final_y/dy)

c = np.zeros((2,xn,yn)) # sets the initial conditions

# adding the glucose on the side 



maximum = 21 # setting an arbitrary maximum because the data is already all 0 at this time
loop_counter = 1
while maximum > 20:
    
    c[1] = c[0]   
    
    for i in range(xn):
        for j in range(yn):
            # corners
            if i == 0 and j == 0: # bottom left
                 c[1, i, j] = 2*Fo_m*(c[0,0,1] + c[0,1,0] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,0,0]
            elif i == 0 and j == yn-1: # upper left
                c[1, i, j] = 2*Fo_m*(c[0,0,(yn-1)-1] + c[0,1,(yn-1)] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,0,(yn-1)]
            elif i == xn-1 and j == 0: # bottom right
                c[1, i, j] = 2*Fo_m*(c[0,(xn-1)-1,0] + c[0,(xn-1),1] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,(xn-1),0]
            elif i == xn-1 and j == yn-1: # upper right
                c[1, i, j] = 2*Fo_m*(c[0,(xn-1),(yn-1)-1] + c[0,(xn-1)-1,(yn-1)] + 2*Sh*C_e) + (1 - 4*Fo_m - 4*Fo_m*Sh)*c[0,(xn-1),(yn-1)]
            
            # sides
            elif i == 0 and (1 <= j <= (yn-1)-1): # left side (exluding corners)
                c[1, i, j] = 2*Fo_m*(c[0,i+1,j] + 1/2*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,i,j]
            
            elif i == xn-1 and (1 <= j <= (yn-1)-1): # right side (exluding corners)
                c[1, i, j] = 2*Fo_m*(c[0,i-1,j] + 1/2*(c[0,i,j-1] + c[0,i,j+1]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,i,j]
            
            elif j == 0 and (1 <= i <= (xn-1)-1): # bottom side (exluding corners)
                c[1, i, j] = 2*Fo_m*(c[0,i,j+1] + 1/2*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,i,j]
            
            elif j == yn-1 and (1 <= i <= (xn-1)-1): # top side (exluding corners)
                c[1, i, j] = 2*Fo_m*(c[0,i,j-1] + 1/2*(c[0,i-1,j] + c[0,i+1,j]) + Sh*C_e) + (1 - 4*Fo_m - 2*Fo_m*Sh)*c[0,i,j]
            
            # interior
            else:
                
                if loop_counter == 1:
                    c[1,i,j] = C_e
                    print(c[1])
                else:
                    c[1, i, j] = c[0,i,j] + D*dt*((c[0,i-1,j] - 2*c[0,i,j] + c[0,i+1,j])/(dx**2) + (c[0,i,j-1] - 2*c[0,i,j] + c[0,i,j+1])/(dy**2))
    
    loop_counter += 1
    maximum = max(max(row) for row in c[1])


print(c[1])
# c) Console Output






