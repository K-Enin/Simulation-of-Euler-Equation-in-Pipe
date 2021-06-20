"""
MAIN CODE
Created on Sat Jun 12 19:09:06 2021

@author: katharinaenin
Isothermal Euler and Weymouth Equation with Lax Friedrich Scheme
Initial Conditions: Lemma 1 in Paper "The isothermal Euler equation for ideal gas
with source term"

The purpose of this code is to compare the LxF solution to the graphics in
the paper.
"""
# Gleichungen mit Quellterm
# Simuliere 1.6 Sekunden

import numpy as np
import matplotlib.pyplot as plt
from math import exp

a_square = 1
teta = 0.7
CFL = 0.8
dx = 0.0625 #40 space steps
n = 40 #n = 2.5/0.0625
T = 2
X = 2.5
beta = -1

# Euler Gleichung
def EulerGleichung():
    # initial conditions
    P_Euler = np.zeros((n,1)) #Spaltenvektor
    P_Euler_new = np.zeros((n,1)) #Spaltenvektor
    
    x = 0
    for i in range(0,n):
        P_Euler[i,0] = exp(beta*x) #Spaltenvektor
        x = x + dx

    Q_Euler = np.zeros((n,1))  
    Q_Euler_new = np.zeros((n,1))      
            
    T = 0
    GeneralStep = 0
    
    while T < 2:
        
        # Q_Euler/P_Euler computes the velocity
        # Calculate dt with the CFL condition
        dt = (CFL*dx)/(np.amax(Q_Euler/P_Euler + a_square*np.ones((1,n))))
        
        for j in range(1,n-1):
            P_Euler_new[j,0] = 0.5*(P_Euler[j-1,0]+P_Euler[j+1,0]) - dt/(2*dx)*(Q_Euler[j+1,0]-Q_Euler[j-1,0])
            Q_Euler_new[j,0] = 0.5*(Q_Euler[j-1,0]+Q_Euler[j+1,0]) - dt/(2*dx)*((Q_Euler[j+1,0]*Q_Euler[j+1,0])/(P_Euler[j+1,0]) - a_square*P_Euler[j+1,0] - (Q_Euler[j-1,0]*Q_Euler[j-1,0])/(P_Euler[j-1,0])+a_square*P_Euler[j-1,0]) - dt*teta/2*(Q_Euler[j-1,0]*abs(Q_Euler[j-1,0])/P_Euler[j-1,0] + Q_Euler[j+1,0]*abs(Q_Euler[j+1,0])/P_Euler[j+1,0]) 
            #0.5*(Q[j-1]+Q[j+1]) - dt/(2*dx)*a_square*(P[j+1,t]-P[j-1,t]) - dt*Lambda/(4*D)*(Q[j-1,t]*abs(Q[j-1,t])/P[j-1,t] + Q[j+1,t]*abs(Q[j+1,t])/P[j+1,t])
        
        
        #Kompatibilitätsbedingungen
        P_Euler_new[j,0] = exp(beta*0)
        P_Euler_new[n-1,0] = exp(beta*2.5)

        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        # print every 300 general steps
        T = T + dt
        GeneralStep = GeneralStep + 1
        
        if GeneralStep%50 == 0:
            print("General Step:" + str(GeneralStep) + " ,Time: " + str(T))
            
    return GeneralStep, dt, T, P_Euler, Q_Euler


if __name__ == '__main__':

    steps1, dt_Euler, T, P_Iso_Euler, Q_Iso_Euler = EulerGleichung()
    P_Iso_Euler_T = P_Iso_Euler.reshape(-1,1)
    Q_Iso_Euler_T = Q_Iso_Euler.reshape(-1,1)
    
    plt.figure(1)
    plt.plot(np.arange(0,n),P_Iso_Euler_T)
    plt.ylabel('Density [bar]')
    plt.xlabel('x_steps')
    plt.title('Simulation of Density for Euler with 0.16 sec')
    
    plt.figure(3)
    plot2 = plt.plot(np.arange(0,n),Q_Iso_Euler_T)
    plt.ylabel('Flux Q')
    plt.xlabel('x_steps')
    plt.title('Simulation of Flux for Euler with 0.16 sec')
    
    plt.show() 