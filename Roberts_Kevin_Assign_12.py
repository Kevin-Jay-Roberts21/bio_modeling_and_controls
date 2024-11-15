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
C = B.T.flatten()
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

A0 = np.zeros((5, 5))
B0 = np.zeros((1, 5))
B0 = B0.T.flatten()

# populating teh A0 and B0 matrices
for i in range(5):
    
    # define terms in the A0 matrix
    if i == 0:
        pass
    elif i == 4:
        pass
    else:
        pass
    
    # define the terms in the B0 matrix
    
    
    
    pass
