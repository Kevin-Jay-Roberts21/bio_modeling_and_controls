#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:32:29 2024

@author: kevinjayroberts
"""

import matplotlib.pyplot as plt
import numpy as np

#############
# PROBELM 1 #
#############

# a)
A = np.arange(10)
B = A.reshape(5,2)
C = B.T
D = A.reshape(2,5)
E = np.array([[2, 1], [5, 4]])
F = np.array([1, 2, 3, 4]).reshape(2, 2)
G = np.dot(E, F)
X = np.linalg.solve(E, G)
I = np.array([[1, 2, 3, 4], [2, 4, 5, 6], [7, 8, 5, -2], [4, 1, 2, 3]])
J = np.array([[5], [8], [1], [5]])
X2 = np.linalg.solve(I, J)

print("PROBLEM 1 OUTPUT:")
print("A: " + str(A))
print()
print("B: " + str(B))
print()
print("C: " + str(C))
print()
print("D: " + str(D))
print()
print("E: " + str(E))
print()
print("F: " + str(F))
print()
print("G: " + str(G))
print()
print("X: " + str(X))
print()
print("I: " + str(I))
print()
print("J: " + str(J))
print()
print("X: " + str(X2))
print()

r_list = [0.1, 0.25, 0.5]
A_list = []
B_list = []
T = [0, 0.2, 0.4, 0.2, 0]


for r in r_list:
    
    A0 = np.zeros((5, 5))
    B0 = np.zeros((5, 1))
    
    # populating teh A0 and B0 matrices
    for i in range(5):
        
        # define terms in the A0 matrix
        if i == 0: # for the first row
            A0[0, 0] = 2 + 2*r
            A0[0, 1] = -r
            B0[0] = r*T[0] + (2 - 2*r)*T[1] + r*T[2] + r*T[0]
        elif i == 4: # for the last row
            A0[4, 4] = 2 + 2*r
            A0[4, 4-1] = -r
            B0[4] = r*T[2] + (2 - 2*r)*T[3] + r*T[4] + r*T[4]
        else: # for the middle rows
            A0[i, i-1] = -r    
            A0[i, i] = 2 + 2*r
            A0[i, i+1] = -r
            B0[i] = r*T[i-1] + (2 - 2*r)*T[i] + r*T[i+1]

    A_list.append(A0)
    B_list.append(B0)

for i in range(len(A_list)):
    print("A0 (r=" + str(r_list[i]) + "): ")
    print(A_list[i])
    
for i in range(len(B_list)):
    print("B0 (r=" + str(r_list[i]) + "): ")
    print(B_list[i])
print()    

X = np.linalg.solve(A_list[2], B_list[2])
print("X (r=0.5): " + str(X))
print()
X_r = np.flip(X, axis=0)
print("X_r: " + str(X_r))
print()
print("X_T: " + str(np.delete(X, 0)))
print()
# js = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# print(np.concatenate((js[0], X, X_r, js[10])))

#############
# PROBLEM 2 #
#############

# defining the parameters
alpha = 1
L = 1
dx = 0.1
dt_max = dx**2/2
dt = dt_max
r = alpha*dt/dx**2

# picking a total time
total_time = 20 # (seconds)

xn = int(L/dx) + 1
tn = int(total_time/dt) + 1

temp = np.zeros((tn, xn))

# setting the initial temp
for i in range(xn):
    if 0 < i*dx <= 0.5:
        temp[0, i] = 2*i*dx
    elif 0.5 <= i*dx < 1:
        temp[0, i] = 2*(1-i*dx)

# running the for loop
for i in range(1, tn): # starting at 1 because we already have the initial condition
    
    # define the A and B matrices
    A = np.zeros((xn, xn))
    B = np.zeros((xn, 1))

    for j in range(xn):
        if j == 0: # for the top row
             A[j, j] = 2 + 2*r
             A[j, j+1] = -r
             B[j] = r*temp[i-1, j] + (2 - 2*r)*temp[i-1, j+1] + r*temp[i-1, j+2] + r*temp[i-1, j]         
        
        elif j == xn-1: # for the bottom row
            A[j, j] = 2 + 2*r
            A[j, j-1] = -2*r
            B[j] = r*temp[i-1, j-2] + (2 - 2*r)*temp[i-1, j-1] + r*temp[i-1, j] + r*temp[i-1, j]
        
        else: # for the middle stuff
            A[j, j-1] = -r    
            A[j, j] = 2 + 2*r
            A[j, j+1] = -r
            B[j] = r*temp[i-1, j-1] + (2 - 2*r)*temp[i-1, j] + r*temp[i-1,j+1]
        
    # now solving for the next time step
    temp[i, :] = np.linalg.solve(A, B).T
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    