# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:52:29 2024

@author: Kevin Roberts
"""

import matplotlib.pyplot as plt
import numpy as np

#############
# PROBLEM 1 #
#############

# part a)
B0 = [[2], [6], [0]]
A0 = [[-0, 1, 0],[-0, 0, 2],[-2, 2, 1]]

X = np.linalg.solve(A0, B0)

print("PROBLEM 1a OUTPUT")
print("A0: " + str(A0))
print()
print("B0: " + str(B0))
print()
print("X: " + str(X))
print()
print()

# part b)
CO = [[2, 6, 0], [0, 0, 2], [1, 0, 2], [0, 2, 1]]
# CO = abs(np.concatenate((B0, A0), axis = 0).T)
AW = [[12.011], [1.008], [15.999]]
MW = np.dot(CO, AW)
MM = X*np.delete(MW,0).reshape(3,1) / MW[0]
print("PROBLEM 1b OUTPUT")
print("CO: " + str(CO))
print()
print("AW: " + str(AW))
print()
print("MW: " + str(MW))
print()
print("MM: " + str(MM))
print()


#############
# PROBLEM 2 #
#############

yeild_of_biomass = 0.48 # (g/g)

CO = np.array([[6, 12, 6, 0], [0, 0, 2, 0], [1, 0, 2, 0], [0, 2, 1, 0], [0, 3, 0, 1], [1, 1.77, 0.49, 0.24], [1, 1.55, 0.31, 0.25]])
AW = np.array([[12.011], [1.008], [15.999], [14.007]])
MW = np.dot(CO, AW)

# HOW TO FIND THESE?
BB0 = 3.45952347
BB1 = 0.78490177

BB = np.array([[BB0], [BB1]])

B0 = np.concatenate((CO[0].reshape(4,1), BB), axis = 0)

cells = np.array([[0, 0, 0, 0, 1, 0]])
prot = np.array([[0, 0, 0, 0, 0, 1]])

A0 = np.concatenate((np.delete(CO, 0, 0).T, cells, prot), axis=0)

X = np.linalg.solve(A0, B0)

MM2 = X*np.delete(MW, 0).reshape(6,1) / MW[0]

print("PROBLEM 2 OUTPUT")
print("MM2: " + str(MM2))



