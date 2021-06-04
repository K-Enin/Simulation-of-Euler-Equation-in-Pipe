#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:42:21 2021
@author: katharinaenin

Isothermal Euler with Simple Upwind Scheme (explicit)
"""

import numpy as np
import math
from numpy import savetxt

# (mxn) Matrix
# m - number of time points
# n - number of space points
m, n = 3600, 75; 

# Define  matrix
P_exp = np.zeros((m,n))
Q_exp = np.zeros((m,n))

# Step sizes & constants
#dt = 1/3600
#dx = 33

dt = 1/60 # (s)
dx = 2000 # (m)

a_square = 300*300
D = 0.5
Lambda = 0.011

# Initial conditions
p_in = 65
q_in = 100

condition = 2

# Test different conditions

# Condition 1 (p,q given at beginning of pipe for all t)
if condition == 1:
    print("Using condition " + str(condition))
    P_exp[0,:] = p_in  #t = 0
    P_exp[:,0] = p_in  
    
    Q_exp[0,:] = q_in  #t = 0
    Q_exp[:,0] = q_in  #at beginning of pipe

# Condition 2 (p given at beginning at q at end of pipe for all t)
elif condition == 2:
    print("Using condition " + str(condition))
    P_exp[0,:] = p_in  #t = 0
    P_exp[:,0] = p_in  #at beginning of the pipe
    
    Q_exp[0,:] = q_in  #t = 0
    Q_exp[:,n-1] = q_in  #at end of the pipe    

# Condition 3 (sudden fall of pressure and flux)
elif condition == 3:
    print("Using condition " + str(condition))
    P_exp[0,:] = p_in  #t = 0
    Q_exp[0,:] = q_in  #t = 0
    
    for i in range(0,m):
        if i <= 300:
            P_exp[i,0] = 60
        else: 
            P_exp[i,0] = 30
    
    for i in range(0,m):
        if i <= 300:
            Q_exp[i,0] = 100
        else: 
            Q_exp[i,0] = 50

# Calculate next steps
if condition == 1 or condition == 3:            
    for t in range(0,m-1):
        for l in range(1,n):
            P_exp[t+1,l] = P_exp[t,l] + dt*((1/dx)*(Q_exp[t,l-1]-Q_exp[t,l]))
            Q_exp[t+1,l] = Q_exp[t,l] - (dt/dx)*(a_square*P_exp[t,l]+Q_exp[t,l]*Q_exp[t,l]/P_exp[t,l]-\
                a_square*P_exp[t,l-1]-Q_exp[t,l-1]*Q_exp[t,l-1]/P_exp[t,l-1])-\
                dt*Lambda*Q_exp[t,l]*abs(Q_exp[t,l])/(2*D*P_exp[t,l])
            
elif condition == 2:
    for t in range(0,m-1):
        for l in range(0,n-1):
            P_exp[t+1,l+1] = P_exp[t,l+1] + (dt/dx)*(Q_exp[t,l]-Q_exp[t,l+1])
        for l in range(1,n):
            #Q_exp[t+1,n-1-l] = Q_exp[t,n-1-l]-(dt/dx)*(a_square*P_exp[t,n-1-l]+Q_exp[t,n-1-l]*Q_exp[t,n-1-l]/P_exp[t,n-1-l]-a_square*P_exp[t,n-l]-Q_exp[t,n-l]*Q_exp[t,n-l]/P_exp[t,n-l])-dt*Lambda*Q_exp[t,n-1-l]*abs(Q_exp[t,n-1-l])/(2*D*P_exp[t,n-1-l])
            Q_exp[t+1,n-l-1] = Q_exp[t,n-l-1]+(dt/dx)*(a_square*P_exp[t,n-l-1]+Q_exp[t,n-l-1]*Q_exp[t,n-l-1]/P_exp[t,n-l-1]-a_square*P_exp[t,n-l]-Q_exp[t,n-l]*Q_exp[t,n-l]/P_exp[t,n-l])-dt*(Lambda*Q_exp[t,n-l-1]*abs(Q_exp[t,n-l-1]))/(2*D*P_exp[t,n-l-2])
        
    savetxt('IsothermalEuler_matrixP.csv', P_exp, fmt = '%10.4f', delimiter = ';')
    savetxt('IsothermalEuler_matrixQ.csv', Q_exp, fmt = '%10.4f', delimiter = ';')