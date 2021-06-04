#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 12:11:08 2021
@author: katharinaenin

Isothermal Euler with implicit Box-Scheme
"""

import numpy as np
import math
from sympy import symbols, Eq, solve

# (mxn) Matrix
# m - number of time points
# n - number of space points
m, n = 60, 75; 

# Define  matrix
P_impl = np.zeros((m,n))
Q_impl = np.zeros((m,n))

# Step sizes & constants
dt = 15 # old: 1s
dx = 2 # 2000(m) klappt nicht

a_square = 115600
D = 0.5
Lambda = 0.011

condition = 3

# Initial conditions
p_in = 60 
q_in = 100

# Test different conditions

# Condition 1
if condition == 1:
    P_impl[0,:] = p_in  #t = 0
    P_impl[:,0] = p_in  #at beginning of pipe
    
    Q_impl[0,:] = q_in  #t = 0
    Q_impl[:,0] = q_in  #at beginning of pipe

# Condition 2 (sudden fall of pressure and flux)
elif condition == 2: 
    P_impl[0,:] = p_in  #t = 0
    Q_impl[0,:] = q_in  #t = 0
    
    for i in range(0,m):
        if i <= 300:
            P_impl[i,0] = 60
        else: 
            P_impl[i,0] = 30
    
    for i in range(0,m):
        if i <= 300:
            Q_impl[i,0] = 100
        else: 
            Q_impl[i,0] = 50

# Condition 3 (smooth fall of pressure)
elif condition == 3:
    P_impl[0,:] = p_in  #t = 0
    Q_impl[0,:] = q_in  #t = 0
    Q_impl[:,0] = q_in
    
    # at beginning of pipe pressure is diminishing
    for i in range(0,m):
        P_impl[i,0] = p_in*math.exp(1/1000*(-i))

# Solve Equation
x, y = symbols('x y')

Dictionary = {}

for t in range(0,m-1):
    for l in range(1,n):
        eq1 = Eq((P_impl[t+1,l-1]+x)/2-(P_impl[t,l-1]+P_impl[t,l])/2+(dt/dx)*(y-Q_impl[t+1,l-1]))
        eq2 = Eq((Q_impl[t+1,l-1]+y)/2-(Q_impl[t,l-1]+Q_impl[t,l])/2+
                 dt/dx*(y*y/x+a_square*x-Q_impl[t+1,l-1]*Q_impl[t+1,l-1]/P_impl[t+1,l-1]-
                        a_square*P_impl[t+1,l-1])-dt/2*(-Lambda*y*y/(2*D*x)-
                            Lambda*Q_impl[t+1,l-1]*Q_impl[t+1,l-1]/(2*D*P_impl[t+1,l-1])))
        sol_dict = solve((eq1,eq2), (x, y)) # here multiple solutions
        
        # find the p, q that are more close to the previos p, q, otherwise overflow
        if abs(P_impl[t,l] - sol_dict[0][0]) > abs(P_impl[t,l] - sol_dict[1][0]):
            P_impl[t+1,l] = sol_dict[1][0]
            Q_impl[t+1,l] = sol_dict[1][1]
        else:
            P_impl[t+1,l] = sol_dict[0][0]
            Q_impl[t+1,l] = sol_dict[0][1]
                
        Dictionary[t+l] = sol_dict
        
        print(sol_dict)
        print("This is P: " + str(P_impl[t+1,l]))
        print("This is Q: " +  str(Q_impl[t+1,l]))