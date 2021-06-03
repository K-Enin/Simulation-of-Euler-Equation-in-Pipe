#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 13:34:37 2021
@author: katharinaenin
"""
import numpy as np
import math

# Weymouth Equation (simplified Isothermal) with Simple Upwind Scheme (explicit)

# (mxn) Matrix
# m - number of time points
# n - number of space points
m, n = 3600, 75; 

# Define  matrix
P = np.zeros((m,n))
Q = np.zeros((m,n))

# Step sizes & constants
# dt = 1/3600
# dx = 33

dt = 1/60 # (s)
dx = 2000 # (m)

a_square = 115600
D = 0.5
Lambda = 0.011

# Initial conditions
p_in = 60 
q_in = 100

condition = 2

# Test different conditions

# Condition 1
if condition == 1:
    P[0,:] = p_in  #t = 0
    #P[:,0] = p_in  
    
    Q[0,:] = q_in  #t = 0
    #Q[:,0] = q_in  #at beginning of pipe

# Condition 2
elif condition == 2:
    P[0,:] = p_in  #t = 0
    P[:,0] = p_in  #at beginning of the pipe
    
    Q[0,:] = q_in  #t = 0
    Q[:,n-1] = q_in  #at end of the pipe
  
# Condition 3
elif condition == 3: 
    
    P[0,:] = p_in  #t = 0
    Q[0,:] = q_in  #t = 0
    
    for i in range(0,m):
        if i <= 300:
            P[i,0] = 60
        else: 
            P[i,0] = 30
    
    for i in range(0,m):
        if i <= 300:
            Q[i,0] = 100
        else: 
            Q[i,0] = 50
    # at beginning of pipe pressure is diminishing
    # for i in range(0,m):
    #    P[i,0] = p_in*math.exp(1/1000*(-i))


# Calculate the values

if condition == 1 or condition == 3:
    for t in range(0,m-1):
        for l in range(1,n-1):
            P[t+1,l] = P[t,l] + dt*((1/dx)*(Q[t,l-1]-Q[t,l]))
            Q[t+1,l] = Q[t,l] + dt*((a_square/dx)*(P[t,l-1]-P[t,l])-Lambda*Q[t,l]*abs(Q[t,l])/(2*D*P[t,l]))


if condition == 2:
    for t in range(0,m-1):
        for l in range(0,n-1):
            P[t+1,l+1] = P[t,l+1] + dt*((1/dx)*(Q[t,l]-Q[t,l+1]))
        for l in range(1,n):
            Q[t+1,n-1-l] = Q[t,n-1-l] + dt*((a_square/dx)*(P[t,n-1-l]+P[t,n-1-l+1])-Lambda*Q[t,n-1-l]*abs(Q[t,n-1-l])/(2*D*P[t,n-1-l]))
