#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:42:21 2021

@author: katharinaenin
"""

import numpy as np
import math

# Isothermal Euler with Simple Upwind Scheme (explicit)

# (mxn) Matrix
# m - number of time points
# n - number of space points
m, n = 3600, 100; 

# Define  matrix
P_exp = np.zeros((m,n))
Q_exp = np.zeros((m,n))

# Step sizes & constants
#dt = 1/3600
#dx = 33

dt = 1/60
dx = 2000

a_square = 115600
D = 0.5
Lambda = 0.011

# Initial conditions
p_in = 60 
q_in = 100

condition = 1

# Test different conditions
if condition == 1:
    ###################
    ### Condition 1 ###
    ###################
    P_exp[0,:] = p_in  #t = 0
    P_exp[:,0] = p_in  
    
    Q_exp[0,:] = q_in  #t = 0
    Q_exp[:,0] = q_in  #at beginning of pipe

elif condition == 2: 
    ###################
    ### Condition 2 ###
    ###################
    
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
    # at beginning of pipe pressure is diminishing
    # for i in range(0,m):
    #    P[i,0] = p_in*math.exp(1/1000*(-i))

# Calculate next steps
if condition == 1 or condition == 2:            
    for t in range(0,m-1):
        for l in range(1,n):
            P_exp[t+1,l] = P_exp[t,l] + dt*((1/dx)*(Q_exp[t,l-1]-Q_exp[t,l]))
            Q_exp[t+1,l] = Q_exp[t,l] - (dt/dx)*(a_square*P_exp[t,l]+Q_exp[t,l]*Q_exp[t,l]/P_exp[t,l]
                               - a_square*P_exp[t,l-1]-Q_exp[t,l-1]*Q_exp[t,l-1]/P_exp[t,l-1])
            -dt*Lambda*Q_exp[t,l]*abs(Q_exp[t,l])/(2*D*P_exp[t,l])