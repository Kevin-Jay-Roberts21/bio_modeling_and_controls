#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:36:01 2024

@author: kevinjayroberts
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# defining the global parameters

T_f = 20 # (C) # target temp
T_i = 4 # (C) # initial temp 
save_step = 10 # (save every 10 seconds)

# other constant params
p = 998.2 # (kg/m^3) desnity
C_p = 4182 # (J/(kg*K)) heat capacity
kh = 0.598 # (W/(m*k)) thermal conductivity
h = 750 # (W/(K*m^2)) heat transfer coefficient
Te = 90 # (C) external temperature
alpha = kh/(p*C_p)

dx = 0.001 # (m)
dy = 0.001 # (m)
dz = 0.001 # (m)
dt = 0.1 # (s)

def run_heat_transfer(x, y, z):
    
    interval_steps = int(save_step / dt)
    
    # defining the lengths of the width, height and thickness, respectively
    xn = int(x/dx)
    yn = int(y/dy)
    zn = int(z/dz)
    
    # calculating the biot number and thermal diffusivity
    Bi = h*dx/kh # biot number
    Td = alpha # thermal diffusivity
    
    # defining the boundary conditions
    T_space_0 = 0 # (C)

    # defining the cube (which will change at every time step)
    c = np.zeros((2, zn+1, yn+1, xn+1))
    
    # zero is already every in c at this step, so we don't need to use T_space_0
    # but now we need to at 4 everywhere in the cube except for the boundary
    c[0,:,:,:] = T_i 
    
    # defining the data (which will be added to at every 10 seconds)
    # it is set at a length of "1" for now, but we will add to it every 10 seconds
    # to add another element to it, we just grab c[1], and do data.append(c[1])
    data = np.zeros((1, zn+1, yn+1, xn+1))
    data[0] = c[0]
    
    step_saved = 0
    loop_counter = 0
    save_step_counter = 0
    minimum = 0 # just needs to be less that T_f to begin with, it will change in the loop
    time = 0
    
    while minimum < T_f:
        
        c[0,:,:,:] = c[1,:,:,:]   
        
        for i in range(0,xn+1): # width loop (x)
            for j in range(0,yn+1): # height loop (y)
                for k in range(0,zn+1): # thickness loop (z)
                    
                
                    # defining cases all of the corners    
                    # left, top, forward
                    if i == 0 and j == 0 and k == 0: 
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    
                    elif i == 0 and j == 0 and k == zn: # left, top, back
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif i == xn and j == 0 and k == zn: # right, top, back
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif i == xn and j == 0 and k == 0: # right, top, forward
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    ##############################################                    
                    elif i == 0 and j == yn and k == 0: # left, bottom, forward
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif i == 0 and j == yn and k == zn: # left, bottom, back
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif i == xn and j == yn and k == zn: # right, bottom, back
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif i == xn and j == yn and k == 0: # right, bottom, forward
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    # EDGES                    
                    elif i == 0 and k == 0: # left, forward edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif i == 0 and k == zn: # left, back edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif i == xn and k == zn: # right, back edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif i == xn and k == 0: # right, forward edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i] )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    #####################################
                    elif j == 0 and k == 0: # front, upper edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif j == 0 and k == zn: # back, upper edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif j == yn and k == zn: # back, lower edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif j == yn and k == 0: # front, lower edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0 )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (1*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 1*c[0,k,j,i-1] )/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))
                    ########################################                    
                    elif i == 0 and j == 0: # left, upper edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif i == xn and j == 0: # right, upper edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif i == xn and j == yn: # right, lower edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif i == 0 and j == yn: # left, lower edge
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (1*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 1*c[0,k-1,j,i] )/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0 )/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0 )/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    # FACES                    
                    elif k == 0: # front face
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0)/ dz**2 +
                            (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                            (c[0,k,j,i+1] - 2 * c[0,k,j,i] + c[0,k,j,i-1])/ dx**2) +
                                2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    elif k == zn: # back face
                         c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                             (2*c[0,k-1,j,i] - 2 * c[0,k,j,i] + 0)/ dz**2 +
                             (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                             (c[0,k,j,i+1] - 2 * c[0,k,j,i] + c[0,k,j,i-1])/ dx**2) +
                                 2*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                 0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                 0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                                      
                    ###################################                    
                    elif i == 0: # left face
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (c[0,k+1,j,i] - 2 * c[0,k,j,i] + c[0,k-1,j,i])/ dz**2 +
                            (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                            (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0)/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif i == xn: # right face
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                            (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                            (2*c[0,k,j,i-1] - 2 * c[0,k,j,i] + 0)/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                    
                    ########################################                    
                    elif j == 0: # top face
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                            (2*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 0)/ dy**2 +
                            (c[0,k,j,i+1] - 2 * c[0,k,j,i] +c[0,k,j,i-1])/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                2*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    elif j == yn: # bottom face
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (c[0,k+1,j,i] - 2 * c[0,k,j,i] +c[0,k-1,j,i])/ dz**2 +
                            (2*c[0,k,j-1,i] - 2 * c[0,k,j,i] + 0)/ dy**2 +
                            (c[0,k,j,i+1] - 2 * c[0,k,j,i] +c[0,k,j,i-1])/ dx**2) +
                                0*dt * h* alpha/(kh*dz)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dy)*(Te-c[0,k,j,i] ) +
                                0*dt * h* alpha/(kh*dx)*(Te-c[0,k,j,i] ))                        
                    # CENTER                            
                    else:
                        c[1,k,j,i] = c[0,k,j,i] + (dt * alpha * (
                            (c[0,k+1,j,i] - 2 * c[0,k,j,i] + c[0,k-1,j,i])/ dz**2 +
                            (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/ dy**2 +
                            (c[0,k,j,i+1] - 2 * c[0,k,j,i] + c[0,k,j,i-1])/ dx**2) )
                    
        loop_counter += 1
        time += dt
        
        # saving the data every 15 steps
        if loop_counter % interval_steps == 0:
            step_saved += 1
            data = np.append(data, np.zeros((1, zn+1, yn+1, xn+1)), axis=0)
            data[step_saved] = c[1,:,:,:]
            minimum = np.min(c[1,:,:,:])
            # print(f"Day {int{time}}: Minimum mass concentration = {minimum:.4f} g/L")
            

        # Check if all agar has at least target concentration
        if np.all(c[1,:,:] >= T_f):
            # print(f"All mass got to at least {C_f} g/L at day {time:.2f}")
            step_saved += 1
            data = np.append(data, np.zeros((1, zn+1, yn+1, xn+1)), axis=0)
            data[step_saved] = c[1,:,:,:]
            end_time = dt*loop_counter
            break
    
    end_time = round(dt*loop_counter,3)
    
    all_data = [data, c, end_time, Bi, Td, alpha]
    
    return all_data

def plot_min_heats(time_grid, dimless_time_grid, mins, cube_id):
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    # First subplot for the original time plot
    ax1.plot(time_grid, mins, marker='o', label='Cube ' + cube_id + ' Heat minimums')

    # Customizing the first subplot
    ax1.set_title('Heat Cube ' + cube_id + ' Mins at every 10 seconds')
    ax1.set_xlabel("Time (seconds)")
    ax1.set_ylabel('Minimum Heat in Cube ' + cube_id + ' Values')
    ax1.set_xticks(time_grid)  # Set x-ticks to match the time values
    ax1.legend()
    ax1.grid()

    # Second subplot for the dimensionless time plot
    ax2.plot(dimless_time_grid, mins, marker='o', label='Cube ' + cube_id + ' Heat minimums')

    # Customizing the second subplot
    ax2.set_title('Heat Cube ' + cube_id + ' Mins (Dimensionless Time)')
    ax2.set_xlabel("Dimensionless Time")
    ax2.set_ylabel('Minimum Heat in Cube ' + cube_id + ' Values')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_xticks(dimless_time_grid)  # Set x-ticks to match the time values
    ax2.legend()
    ax2.grid()

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

def display_heat_data(cube_data, cube_lengths, cube_id):
    print("#################")
    print("Results of Cube " + cube_id)
    print("x: " + cube_lengths + " cm")
    print("y: " + cube_lengths + " cm")
    print("z: " + cube_lengths + " cm")
    print("Thermal diffusivity: " + str(cube_data[4])) 
    print("Biot Number: " + str(round(cube_data[3], 3)))
    
    # printing the minimum concentration at a certain position for every 10 seconds
    for i in range(len(cube_data[0])):
        
        if i == len(cube_data[0])-1:
            time_of_min_heat = str(cube_data[2])
        else:
            time_of_min_heat = str(i*10)
        
        position = np.unravel_index(np.argmin(cube_data[0][i]), cube_data[0][i].shape)
        print("T = " + time_of_min_heat + " seconds:   min heat: " + str(round(np.min(cube_data[0][i]), 3)) + "    position: " + str(position))


# defining the cube 1, 2, and 3 dimensions
cube1_width = 0.01 # (m) the x value
cube1_height = 0.01 # (m) the y value
cube1_thickness = 0.01 # (m) the z value

cube2_width = 0.011 # (m) the x value
cube2_height = 0.011 # (m) the y value
cube2_thickness = 0.011 # (m) the z value

cube3_width = 0.012 # (m) the x value
cube3_height = 0.012 # (m) the y value
cube3_thickness = 0.012 # (m) the z value

cube1_data = run_heat_transfer(cube1_width, cube1_height, cube1_thickness)
cube2_data = run_heat_transfer(cube2_width, cube2_height, cube2_thickness)
cube3_data = run_heat_transfer(cube3_width, cube3_height, cube3_thickness)

print("PROBLEM 1a OUTPUT:")
display_heat_data(cube1_data, "1", "1")
display_heat_data(cube2_data, "1.1", "2")
display_heat_data(cube3_data, "1.2", "3")

print()
print()


# b) Figure 1a

# getting the all the minimum values at each 10 days for each cube
cube1_mins = []
cube2_mins = []
cube3_mins = []

for i in range(len(cube1_data[0])):
    cube1_mins.append(np.min(cube1_data[0][i]))
for i in range(len(cube2_data[0])):
    cube2_mins.append(np.min(cube2_data[0][i]))
for i in range(len(cube3_data[0])):
    cube3_mins.append(np.min(cube3_data[0][i]))

cube1_L = cube1_width # (m)
cube2_L = cube2_width # (m)
cube3_L = cube3_width # (m)

# this time grids may be longer or shorter for each cube
cube1_t = np.arange(0, (len(cube1_data[0])-1)*save_step, save_step)
cube1_t = np.append(cube1_t, cube1_data[2])

cube2_t = np.arange(0, (len(cube2_data[0])-1)*save_step, save_step)
cube2_t = np.append(cube2_t, cube2_data[2])

cube3_t = np.arange(0, (len(cube3_data[0])-1)*save_step, save_step)
cube3_t = np.append(cube3_t, cube3_data[2])

cube1_dim_less_time = [round(i * (cube1_data[5] / cube1_L**2), 3) for i in cube1_t]
cube2_dim_less_time = [round(i * (cube2_data[5] / cube2_L**2), 3) for i in cube2_t]
cube3_dim_less_time = [round(i * (cube3_data[5] / cube3_L**2), 3) for i in cube3_t]

plot_min_heats(cube1_t, cube1_dim_less_time, cube1_mins, "1") # plotting cube 1 mins
plot_min_heats(cube2_t, cube2_dim_less_time, cube2_mins, "2") # plotting cube 2 mins
plot_min_heats(cube3_t, cube3_dim_less_time, cube3_mins, "3") # plotting cube 3 mins

    
# c) Explain
print("PROBLEM 1b OUTPUT:")
print("i) The differences between panel 1 for each cube, we see that as the heat increases, it takes longer for the heat to get to 20 C.")
print("ii) The panels between i and ii are very similar for each cube, only the time scale is different.")
print("iii) The differences between panel ii are the same types of differences we see in the panel i differences, i.e. it just takes longer for the heat to get to 20 C for the larger cubes.")
print("iv) The implications are that it takes longer for the heat to migrate through the cube if it is larger.")
print()
print()


# d) Figure 2b
##############
# FOR CUBE 1 #
##############
#generate the 'n' plots from the 'n' planes of data in the storage array
xn = int(cube1_width/dx) + 1
yn = int(cube1_height/dy) + 1
for k in range(len(cube1_data[0])):  # Iterate through each cube
    fig = plt.figure()
    
    # Initialize x-axis and y-axis for the plot with correct 2D shape
    xaxis_1d = np.linspace(0, (xn - 1) * dx * 100, xn)  # Create x-axis in centimeters
    yaxis_1d = np.linspace(0, (yn - 1) * dy * 100, yn)  # Create y-axis in centimeters

    xaxis, yaxis = np.meshgrid(xaxis_1d, yaxis_1d)

    # Plotting
    ax = plt.axes(projection='3d')
    ax.set_box_aspect(aspect=(1, 2, 1))  # Scale axes to match ranges

    # Specify z level to slice (e.g., 0, 1, or 2)
    z_level = 1  # Use any fixed z level you need (0, 1, or 2)
    z = cube1_data[0][k, z_level, :, :]  # Select only the 2D slice for the chosen z level

    # Plot surface with x-y plane data
    surf = ax.plot_surface(xaxis, yaxis, z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
    
    fig.colorbar(surf)  # Insert a color bar on the graph
    
    if i == len(cube1_data[0])-1:
        heat_time = cube1_data[2]
    else:
        heat_time = i*10
    
    plt.title(f"Bottom side of Cube 1 at time {heat_time}, z-level {z_level}")  # Title includes cube index and z-level
    plt.xlabel("X-axis (cm)")
    plt.ylabel("Y-axis (cm)")
    
    plt.show()

##############
# FOR CUBE 2 #
##############

xn = int(cube2_width/dx) + 1
yn = int(cube2_height/dy) + 1
for k in range(len(cube2_data[0])):  # Iterate through each cube
    fig = plt.figure()
    
    # Initialize x-axis and y-axis for the plot with correct 2D shape
    xaxis_1d = np.linspace(0, (xn - 1) * dx * 100, xn)  # Create x-axis in centimeters
    yaxis_1d = np.linspace(0, (yn - 1) * dy * 100, yn)  # Create y-axis in centimeters

    xaxis, yaxis = np.meshgrid(xaxis_1d, yaxis_1d)

    # Plotting
    ax = plt.axes(projection='3d')
    ax.set_box_aspect(aspect=(1, 2, 1))  # Scale axes to match ranges

    # Specify z level to slice (e.g., 0, 1, or 2)
    z_level = 1  # Use any fixed z level you need (0, 1, or 2)
    z = cube2_data[0][k, z_level, :, :]  # Select only the 2D slice for the chosen z level
    
    # Plot surface with x-y plane data
    surf = ax.plot_surface(xaxis, yaxis, z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
    
    if i == len(cube2_data[0])-1:
        heat_time = cube2_data[2]
    else:
        heat_time = i*10
    
    fig.colorbar(surf)  # Insert a color bar on the graph
    plt.title(f"Bottom side of Cube 2 at time {heat_time}, z-level {z_level}")  # Title includes cube index and z-level
    plt.xlabel("X-axis (cm)")
    plt.ylabel("Y-axis (cm)")
    
    plt.show()

##############
# FOR CUBE 3 #
##############

xn = int(cube3_width/dx) + 1
yn = int(cube3_height/dy) + 1
for k in range(len(cube3_data[0])):  # Iterate through each cube
    fig = plt.figure()
    
    # Initialize x-axis and y-axis for the plot with correct 2D shape
    xaxis_1d = np.linspace(0, (xn - 1) * dx * 100, xn)  # Create x-axis in centimeters
    yaxis_1d = np.linspace(0, (yn - 1) * dy * 100, yn)  # Create y-axis in centimeters

    xaxis, yaxis = np.meshgrid(xaxis_1d, yaxis_1d)

    # Plotting
    ax = plt.axes(projection='3d')
    ax.set_box_aspect(aspect=(1, 2, 1))  # Scale axes to match ranges

    # Specify z level to slice (e.g., 0, 1, or 2)
    z_level = 1  # Use any fixed z level you need (0, 1, or 2)
    z = cube3_data[0][k, z_level, :, :]  # Select only the 2D slice for the chosen z level

    if i == len(cube3_data[0])-1:
        time_of_min_heat = cube3_data[2]
    else:
        time_of_min_heat = i*10

    # Plot surface with x-y plane data
    surf = ax.plot_surface(xaxis, yaxis, z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
    
    fig.colorbar(surf)  # Insert a color bar on the graph
    plt.title(f"Bottom side of Cube 3 at time {heat_time}, z-level {z_level}")  # Title includes cube index and z-level
    plt.xlabel("X-axis (cm)")
    plt.ylabel("Y-axis (cm)")
    
    plt.show() 



# e) Explain
print("PROBLEM 1e OUTPUT:")
print("i) The shape of the T contour in the plot looks somewhat like a pool after the initial heat is change, where the initial heat is constant at 4 C. This is not just the case for the x-y plane but for all the planes.")
print("ii) 3D model is required to monitor T_min in these cubes, because we cannot simply model just one of the sides, there are multiple dependencies.")
print("iii) The cube being insulted on certain surfaces, the constant convection at the bottom of the cube, and temperature in the external space. Also the diffusion.")
print("iv) increased density, increased heat capacity, or increased thermal conductivity.")



