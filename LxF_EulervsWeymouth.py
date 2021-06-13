#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN CODE
Created on Sat Jun 12 19:09:06 2021

@author: katharinaenin
Isothermal Euler and Weymouth Equation with Lax Friedrich Scheme
Initial Conditions: Sud Shock Tube

The purpose of this code is to demonstrate the difference between isothermal 
Euler Equation and its simplified version - Weymouth equation.

"""
# Gleichungen mit Quellterm
# Simuliere 1.6 Sekunden

import numpy as np
import matplotlib.pyplot as plt

Lambda = 0.0011
D = 0.5
a_square = 377.9683*377.9683
CFL = 0.8
tMax = 0.164 #lass max 0.1 Sekunden laufen
dx = 0.125
a = 377.9683
n = 80 #n=10/0.125

# Euler Gleichung
def EulerGleichung():
    # initial conditions
    P_Euler = np.zeros((1,n))
    P_Euler_new = np.zeros((1,n))

    for i in range(0,n):
        if i <= 40:
            P_Euler[0,i] = 1
        else:
            P_Euler[0,i] = 0.125

    Q_Euler = np.zeros((1,n))  
    Q_Euler_new = np.zeros((1,n))      
            
    T = 0
    GeneralStep = 0
    
    while T < tMax:
        
        # Q_Euler/P_Euler computes the velocity
        # Calculate dt with the CFL condition
        dt = (CFL*dx)/(np.amax(Q_Euler/P_Euler + a_square*np.ones((1,n))))
        
        for j in range(1,n-1):
            P_Euler_new[0,j] = 0.5*(P_Euler[0,j-1]+P_Euler[0,j+1]) - dt/(2*dx)*(Q_Euler[0,j+1]-Q_Euler[0,j-1])
            Q_Euler_new[0,j] = 0.5*(Q_Euler[0,j-1]+Q_Euler[0,j+1]) - dt/(2*dx)*((Q_Euler[0,j+1]*Q_Euler[0,j+1])/(P_Euler[0,j+1])+a_square*P_Euler[0,j+1] - (Q_Euler[0,j-1]*Q_Euler[0,j-1])/(P_Euler[0,j-1])-a_square*P_Euler[0,j-1]) - dt*Lambda/(4*D)*(Q_Euler[0,j-1]*abs(Q_Euler[0,j-1])/P_Euler[0,j-1] + Q_Euler[0,j+1]*abs(Q_Euler[0,j+1])/P_Euler[0,j+1]) 
            #0.5*(Q[j-1]+Q[j+1]) - dt/(2*dx)*a_square*(P[j+1,t]-P[j-1,t]) - dt*Lambda/(4*D)*(Q[j-1,t]*abs(Q[j-1,t])/P[j-1,t] + Q[j+1,t]*abs(Q[j+1,t])/P[j+1,t])
        
        #Kompatibilitätsbedingungen
        P_Euler_new[0,0] = 1
        P_Euler_new[0,n-1] = 0.125

        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        # print every 300 general steps
        T = T + dt
        GeneralStep = GeneralStep + 1
        
        if GeneralStep%500 == 0:
            print("General Step:" + str(GeneralStep) + " ,Time: " + str(T))
            
    return GeneralStep, dt, T, P_Euler, Q_Euler
            
# Weymouth Gleichungen
def WeymouthGleichung():

    # initial conditions
    P_Euler = np.zeros((1,n))
    P_Euler_new = np.zeros((1,n))

    for i in range(0,n):
        if i <= 40:
            P_Euler[0,i] = 1
        else:
            P_Euler[0,i] = 0.125

    Q_Euler = np.zeros((1,n))  
    Q_Euler_new = np.zeros((1,n))      
            
    T = 0
    GeneralStep = 0
    dt = 7/10000000 #7*10^(-7)
    
    while T < tMax:
        
        for j in range(1,n-1):
            P_Euler_new[0,j] = 0.5*(P_Euler[0,j-1]+P_Euler[0,j+1])-dt/(2*dx)*(Q_Euler[0,j+1]-Q_Euler[0,j-1]) #gleich
            Q_Euler_new[0,j] = 0.5*(Q_Euler[0,j-1]+Q_Euler[0,j+1])-dt/(2*dx)*a_square*(P_Euler[0,j+1]-P_Euler[0,j-1]) - dt*Lambda/(4*D)*(Q_Euler[0,j-1]*abs(Q_Euler[0,j-1])/P_Euler[0,j-1] + Q_Euler[0,j+1]*abs(Q_Euler[0,j+1])/P_Euler[0,j+1])
        
        #Kompatibilitätsbedingungen
        P_Euler_new[0,0] = 1
        P_Euler_new[0,n-1] = 0.125
        
        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        T = T + dt
        GeneralStep = GeneralStep + 1
        
        if GeneralStep%500 == 0:
            print("General Step:" + str(GeneralStep) + " ,Time: " + str(T))  
            
    return GeneralStep, T, P_Euler, Q_Euler
            
if __name__ == '__main__':

    steps1, dt_Euler, T, P_Iso_Euler, Q_Iso_Euler = EulerGleichung()
    P_Iso_Euler_T = P_Iso_Euler.reshape(-1,1)
    Q_Iso_Euler_T = Q_Iso_Euler.reshape(-1,1)
    
    plt.figure(1)
    plt.plot(np.arange(0,n),P_Iso_Euler_T)
    plt.ylabel('Density [bar]')
    plt.xlabel('x_steps')
    plt.title('Simulation of Density for Euler with 1.6 sec')
    
    plt.figure(2)
    plot2 = plt.plot(np.arange(0,n),Q_Iso_Euler_T/P_Iso_Euler_T)
    plt.ylabel('Velocity')
    plt.xlabel('x_steps')
    plt.title('Simulation of Velocity for Euler with 1.6 sec')
    
    plt.show() # lasse es weg, wenn der untere Abschnitt auch läuft
    
    # steps, T_wey, P_Wey, Q_Wey = WeymouthGleichung()
    # P_Wey_T = P_Wey.reshape(-1,1)
    # Q_Wey_T = Q_Wey.reshape(-1,1)
    
    # plt.figure(1)
    # plt.plot(np.arange(0,n),P_Wey_T)
    # plt.ylabel('Density [bar]')
    # plt.xlabel('x_steps')
    # plt.title('Simulation of Density for Weymouth with 1.6 sec')
    
    # plt.figure(2)
    # plot2 = plt.plot(np.arange(0,n),Q_Wey_T/P_Wey_T)
    # plt.ylabel('Velocity')
    # plt.xlabel('x_steps')
    # plt.title('Simulation of Velocity for Weymouth with 1.6 sec')        
    
    # plt.show()