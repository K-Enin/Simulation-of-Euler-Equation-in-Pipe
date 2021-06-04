#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 20:11:21 2021

@author: katharinaenin
"""
import numpy as np
from numpy import loadtxt
m, n = 3600, 75; 

#data_dict = load_dictionary(FirstResults.spydata)
#with open("Weymouth_matrixP.txt", 'wb') as f:
#    A = np.load(f)
P_exp = loadtxt('IsothermalEuler_matrixP.csv', delimiter = ';')
Q_exp = loadtxt('IsothermalEuler_matrixQ.csv', delimiter = ';')

P_W = loadtxt('Weymouth_matrixP.csv', delimiter = ';')
Q_W = loadtxt('Weymouth_matrixQ.csv', delimiter = ';')

P_diff = P_exp - P_W
Q_diff = Q_exp - Q_W

maximum_deviation_P = abs(P_diff[0,0])
maximum_deviation_Q = abs(Q_diff[0,0])

for i in range(0,m):
    for j in range(0,n):
        if maximum_deviation_P < abs(P_diff[i,j]):
            maximum_deviation_P = abs(P_diff[i,j])
        
        if maximum_deviation_Q < abs(Q_diff[i,j]):
            maximum_deviation_Q = abs(Q_diff[i,j])
            
print("Maximum deviation for P is: " + str(maximum_deviation_P))
print("Maximum deviation for Q is: " + str(maximum_deviation_Q))

# Here: 
# Maximum deviation for P is: 0.0007000000000090267
# Maximum deviation for Q is: 0.17030000000000456