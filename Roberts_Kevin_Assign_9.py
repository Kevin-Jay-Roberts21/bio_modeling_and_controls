#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:28:48 2024

@author: kevinjayroberts
"""

#############
# PROBLEM 1 #
#############
z = 0.4 # length (m)
y = 0.1 # thickness (m)
x = 0.1 # height (m)

D = 0.252*(10**-4) # (m^2/day)
dx = 0.01 # (cm)
dy = 0.01 # (cm)
k_m = 3.1536*10**(-3) 
Sh = k_m*dx/D
interior = dx**2/(4*D)
side = dx**2/(2*D*(2 + Sh))
corner = dx**2/(4*D*(1 + Sh))
dt_maxes = [interior, side, corner] 
dt_max = min(dt_maxes)
dt = round(dt_max, 3)
Fo_m = D*dt/(dx**2)


