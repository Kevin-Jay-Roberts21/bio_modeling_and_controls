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
save_step = 10 # (save every 10 seconds)

# other constant params
p = 998.2 # (kg/m^3) desnity
C_p = 4182 # (J/(kg*K)) heat capacity
kh = 0.598 # (W/(m*k)) thermal conductivity
h = 750 # (W/(K*m^2)) heat transfer coefficient
T_e = 90 # (C) external temperature
alpha = kh/(p*C_p)

dx = 0.001 # (m)
dy = 0.001 # (m)
dz = 0.001 # (m)
dt = 0.1 # (s)

def run_heat_transfer(x, y, z):
    
    maximum = 60 # (maximum time in seconds)
    S = int(maximum / save_step) + 1
    interval_steps = int(save_step / dt)
    
    # defining the lengths of the width, height and thickness, respectively
    xn = int(x/dx) + 1
    yn = int(y/dy) + 1
    zn = int(z/dz) + 1
    
    # calculating the biot number and thermal diffusivity
    Bi = h*dx/kh # biot number
    Td = alpha # thermal diffusivity
    
    # defining the initial conditions
    T_time_0 = 4 # (C)
    
    # defining the boundary conditions
    T_space_0 = 0 # (C)

    # defining the cube (which will change at every time step)
    c = np.zeros((2, zn, yn, xn))
    
    # zero is already every in c at this step, so we don't need to use T_space_0
    # but now we need to at 4 everywhere in the cube except for the boundary
    c[0,:,:,:] = T_time_0 
    
    # defining the data (which will be added to at every 10 seconds)
    # it is set at a length of "1" for now, but we will add to it every 10 seconds
    # to add another element to it, we just grab c[1], and do data.append(c[1])
    data = np.zeros((S, zn, yn, xn))
    data[0] = c[0]
    
    loop_counter = 0
    save_step_counter = 0
    minimum = 0 # just needs to be less that T_f to begin with, it will change in the loop
    time = 0
    
    while time < maximum:
        
        c[0,:,:,:] = c[1,:,:,:]   
        
        for i in range(xn): # width loop (x)
            for j in range(yn): # height loop (y)
                for k in range(zn): # thickness loop (z)
                    
                
                    # defining cases all of the corners    
                
                    # for the upper front left corner
                    if i == 0 and j == 0 and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 + 
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the upper front right corner
                    elif i == xn-1 and j == 0 and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 + 
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the lower front left corner
                    elif i == 0 and j == yn-1 and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 + 
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the lower front right corner
                    elif i == xn-1 and j == yn-1 and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 + 
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the upper back left corner
                    if i == 0 and j == 0 and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 + 
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the upper back right corner
                    elif i == xn-1 and j == 0 and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 + 
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the lower back left corner
                    elif i == 0 and j == yn-1 and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 + 
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the lower back right corner
                    elif i == xn-1 and j == yn-1 and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 + 
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) + 
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) + 
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) + 
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    
                    
                    
                    
                    
                    
                    # defining cases for all the edges
                    
                    #### FRONT ####
                    # for the front left edge
                    elif i == 0 and (0 < j < yn-1) and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2 * c[0,k,j,i] + 0)/dz**2 +
                                          (1*c[0,k,j+1,i] - 2 * c[0,k,j,i] + 1*c[0,k,j-1,i])/dy**2 +
                                          (2*c[0,k,j,i+1] - 2 * c[0,k,j,i] + 0)/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the front right edge
                    elif i == xn-1 and (0 < j < yn-1) and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 +
                                          (1*c[0,k,j+1,i] - 2*c[0,k,j,i] + 1*c[0,k,j-1,i])/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the front top edge
                    elif (0 < i < xn-1) and j == 0 and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 +
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (1*c[0,k,j,i+1] - 2*c[0,k,j,i] + 1*c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the front bottom edge
                    elif (0 < i < xn-1) and j == yn-1 and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (1*c[0,k,j,i+1] - 2*c[0,k,j,i] + 1*c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    
                    
                    #### MIDDLE ####
                    # for the middle top left edge
                    elif i == 0 and j == 0 and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((1*c[0,k+1,j,i] - 2*c[0,k,j,i] + 1*c[0,k-1,j,i])/dz**2 +
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the middle top right edge
                    elif i == xn-1 and j == 0 and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((1*c[0,k+1,j,i] - 2*c[0,k,j,i] + 1*c[0,k-1,j,i])/dz**2 +
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the middle bottom left edge
                    elif i == 0 and j == yn-1 and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((1*c[0,k+1,j,i] - 2*c[0,k,j,i] + 1*c[0,k-1,j,i])/dz**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the middle bottom right edge
                    elif i == xn-1 and j == yn-1 and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((1*c[0,k+1,j,i] - 2*c[0,k,j,i] + 1*c[0,k-1,j,i])/dz**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    
                    
                    #### BACK ####
                    # for the back left edge
                    elif i == 0 and (0 < j < yn-1) and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 +
                                          (1*c[0,k,j+1,i] - 2*c[0,k,j,i] + 1*c[0,k,j-1,i])/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the back right edge
                    elif i == xn-1 and (0 < j < yn-1) and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 +
                                          (1*c[0,k,j+1,i] - 2*c[0,k,j,i] + 1*c[0,k,j-1,i])/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the back top edge
                    elif (0 < i < xn-1) and j == 0 and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 +
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (1*c[0,k,j,i+1] - 2*c[0,k,j,i] + 1*c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the back bottom edge
                    elif (0 < i < xn-1) and j == yn-1 and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (1*c[0,k,j,i+1] - 2*c[0,k,j,i] + 1*c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    
                    
                    
                    
                    # defining the faces
                    
                    # for the front face
                    elif (0 < i < xn-1) and (0 < j < yn-1) and k == 0:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((2*c[0,k+1,j,i] - 2*c[0,k,j,i] + 0)/dz**2 +
                                          (c[0,k,j+1,i] - 2*c[0,k,j,i] + c[0,k,j-1,i])/dy**2 +
                                          (c[0,k,j,i+1] - 2*c[0,k,j,i] + c[0,k,j,i-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the left face
                    elif i == 0 and (0 < j < yn-1) and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((c[0,k+1,j,i] - 2*c[0,k,j,i] + c[0,k-1,j,i])/dz**2 +
                                          (c[0,k,j+1,i] - 2*c[0,k,j,i] + c[0,k,j-1,i])/dy**2 +
                                          (2*c[0,k,j,i+1] - 2*c[0,k,j,i] + 0)/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the right face
                    elif i == xn-1 and (0 < j < yn-1) and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((c[0,k+1,j,i] - 2*c[0,k,j,i] + c[0,k-1,j,i])/dz**2 +
                                          (c[0,k,j+1,i] - 2*c[0,k,j,i] + c[0,k,j-1,i])/dy**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,i,j,k-1])/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the top face
                    elif (0 < i < xn-1) and j == 0 and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((c[0,k+1,j,i] - 2*c[0,k,j,i] + c[0,k-1,j,i])/dz**2 +
                                          (2*c[0,k,j+1,i] - 2*c[0,k,j,i] + 0)/dy**2 +
                                          (c[0,k,j,i+1] - 2*c[0,k,j,i] + c[0,i,j,k-1])/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            2*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the bottom face
                    elif (0 < i < xn-1) and j == yn-1 and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((c[0,k+1,j,i] - 2*c[0,k,j,i] + c[0,k-1,j,i])/dz**2 +
                                          (0 - 2*c[0,k,j,i] + 2*c[0,k,j-1,i])/dy**2 +
                                          (c[0,k,j,i+1] - 2*c[0,k,j,i] + c[0,i,j,k-1])/dx**2) +
                            0*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the back face
                    elif (0 < i < xn-1) and (0 < j < yn-1) and k == zn-1:
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((0 - 2*c[0,k,j,i] + 2*c[0,k-1,j,i])/dz**2 +
                                          (c[0,k,j+1,i] - 2*c[0,k,j,i] + c[0,k,j-1,i])/dy**2 +
                                          (c[0,k,j,i+1] - 2*c[0,k,j,i] + c[0,i,j,k-1])/dx**2) +
                            2*dt*h*alpha/(kh*dz)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dy)*(T_e - c[0,k,j,i]) +
                            0*dt*h*alpha/(kh*dx)*(T_e - c[0,k,j,i]))
                    
                    # for the interior points
                    elif (0 < i < xn-1) and (0 < j < yn-1) and (0 < k < zn-1):
                        c[1,k,j,i] = c[0,k,j,i] + (
                            dt * alpha * ((c[0,k+1,j,i] - 2 * c[0,k,j,i] + c[0,k-1,j,i])/dz**2 +
                                          (c[0,k,j+1,i] - 2 * c[0,k,j,i] + c[0,k,j-1,i])/dy**2 +
                                          (c[0,k,j,i+1] - 2 * c[0,k,j,i] + c[0,k,j,i-1])/dx**2))
                    
                    
                
        # increasing the loop_counter
        loop_counter += 1
        time += dt
        
        # getting the minimum heat of the entire cube at this time step
        minimum = np.min(c[1,:,:,:])
        
        # adding to the data at every 10 seconds
        if loop_counter % interval_steps == 0:
            save_step_counter += 1
            data[save_step_counter] = c[1]
            
            
        # Check if heat has at least target concentration
        if np.all(c[1,:,:,:] >= T_f):
            save_step_counter += 1            
            data[save_step_counter] = c[1]
            end_time = round(dt*loop_counter,3)
            break
            
            
        print("Loop count: " + str(loop_counter))
    
    print("Finished a cube loop")
    end_time = round(dt*loop_counter,3)
    
    all_data = [data, c, end_time, Bi, Td, alpha, save_step_counter+1]
    
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
    for i in range(cube_data[6]):
        
        if i == cube_data[6]-1:
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

for i in range(cube1_data[6]):
    cube1_mins.append(np.min(cube1_data[0][i]))
for i in range(cube2_data[6]):
    cube2_mins.append(np.min(cube2_data[0][i]))
for i in range(cube3_data[6]):
    cube3_mins.append(np.min(cube3_data[0][i]))

cube1_L = cube1_width # (m)
cube2_L = cube2_width # (m)
cube3_L = cube3_width # (m)

# this time grids may be longer or shorter for each cube
cube1_t = [0, 10, 20, 30, int(cube1_data[2])]
cube2_t = [0, 10, 20, 30, int(cube2_data[2])]
cube3_t = [0, 10, 20, 30, 40, int(cube3_data[2])]

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
for k in range(cube1_data[6]):  # Iterate through each cube
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
    
    # Check shapes to ensure compatibility
    print("xaxis shape:", xaxis.shape)  # Should be (3, 3)
    print("yaxis shape:", yaxis.shape)  # Should be (3, 3)
    print("z shape:", z.shape)          # Should be (3, 3)

    # Plot surface with x-y plane data
    surf = ax.plot_surface(xaxis, yaxis, z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
    
    fig.colorbar(surf)  # Insert a color bar on the graph
    
    if i == cube1_data[6]-1:
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
for k in range(cube2_data[6]):  # Iterate through each cube
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
    
    # Check shapes to ensure compatibility
    print("xaxis shape:", xaxis.shape)  # Should be (3, 3)
    print("yaxis shape:", yaxis.shape)  # Should be (3, 3)
    print("z shape:", z.shape)          # Should be (3, 3)

    # Plot surface with x-y plane data
    surf = ax.plot_surface(xaxis, yaxis, z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
    
    if i == cube2_data[6]-1:
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
for k in range(cube3_data[6]):  # Iterate through each cube
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

    if i == cube3_data[6]-1:
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



