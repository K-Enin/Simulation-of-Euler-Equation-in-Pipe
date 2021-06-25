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
from mpl_toolkits import mplot3d

a_square = 1
teta = 0.7
CFL = 0.1 #0.8
dx = 0.0625 #40 space steps
n = 40 #n = 2.5/0.0625
T = 2
X = 2.5
beta = -1


# Euler Gleichung mit LxF
def EulerGleichung_LxF():
    # initial conditions
    P_Euler = np.zeros((n+1,1)) #Spaltenvektor
    P_Euler_new = np.zeros((n+1,1)) #Spaltenvektor
    
    x = 0
    X_array = [x]
    for i in range(0,n+1):
        P_Euler[i,0] = exp(beta*x) #Spaltenvektor
        x = x + dx
        X_array = np.append(X_array,x)
        
    Q_Euler = np.zeros((n+1,1))  #initialisiere
    Q_Euler_new = np.zeros((n+1,1))      
            
    P_Matrix = P_Euler
    Q_Matrix = Q_Euler
    
    T = 0
    GeneralSteps = 0
    T_array = [T]
    
    while T < 2:
        
        # Q_Euler/P_Euler computes the velocity
        # Calculate dt with the CFL condition
        dt = (CFL*dx)/(np.amax(abs(Q_Euler/P_Euler) + a_square*np.ones((1,n+1))))
        #dt = 0.001
        
        # Get the time value 
        T = T + dt
        T_array = np.append(T_array, T)
        
        for j in range(1,n):
            P_Euler_new[j,0] = 0.5*(P_Euler[j-1,0]+P_Euler[j+1,0]) - dt/(2*dx)*(Q_Euler[j+1,0]-Q_Euler[j-1,0])
            Q_Euler_new[j,0] = 0.5*(Q_Euler[j-1,0]+Q_Euler[j+1,0]) - dt/(2*dx)*((Q_Euler[j+1,0]*Q_Euler[j+1,0])/(P_Euler[j+1,0]) + a_square*P_Euler[j+1,0] - (Q_Euler[j-1,0]*Q_Euler[j-1,0])/(P_Euler[j-1,0])-a_square*P_Euler[j-1,0]) - dt*teta/4*(Q_Euler[j-1,0]*abs(Q_Euler[j-1,0])/P_Euler[j-1,0] + Q_Euler[j+1,0]*abs(Q_Euler[j+1,0])/P_Euler[j+1,0]) 
            #0.5*(Q[j-1]+Q[j+1]) - dt/(2*dx)*a_square*(P[j+1,t]-P[j-1,t]) - dt*Lambda/(4*D)*(Q[j-1,t]*abs(Q[j-1,t])/P[j-1,t] + Q[j+1,t]*abs(Q[j+1,t])/P[j+1,t])
        
        #Kompatibilitätsbedingungen
        P_Euler_new[0,0] = 4.75*T*T+1 #2.7067*exp(T) #node 0
        P_Euler_new[n,0] = P_Euler[n,0]         #last node, Kompatibilität

        Q_Euler_new[0,0] = 10*T*T #5.4134*exp(T) #node 0
        Q_Euler_new[n,0] = Q_Euler[n,0]          #last node

        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        P_Matrix = np.append(P_Matrix, P_Euler, axis = 1)
        Q_Matrix = np.append(Q_Matrix, Q_Euler, axis = 1)
        
        # print every 300 general steps
        GeneralSteps = GeneralSteps + 1
        
        if GeneralSteps%50 == 0:
            print("General Step:" + str(GeneralSteps) + " ,Time: " + str(T))
            

    return dt, GeneralSteps, X_array[:-1], T_array, P_Matrix, Q_Matrix


# Weymouth Gleichung mit LxF
def Weymouth_LxF():
    # initial conditions
    P_Euler = np.zeros((n+1,1)) #Spaltenvektor
    P_Euler_new = np.zeros((n+1,1)) #Spaltenvektor
    
    x = 0
    X_array = [x]
    for i in range(0,n+1):
        P_Euler[i,0] = exp(beta*x) #Spaltenvektor
        x = x + dx
        X_array = np.append(X_array,x)
        
    Q_Euler = np.zeros((n+1,1))  #initialisiere
    Q_Euler_new = np.zeros((n+1,1))      
            
    P_Matrix = P_Euler
    Q_Matrix = Q_Euler
    
    T = 0
    GeneralSteps = 0
    T_array = [T]
    
    # in order to fulfill CFL condition, we need 
    # dt*a <= dx -> dt*1 <= 0.0625 -> dt = 0.05
    dt = 0.001
    
    while T < 2:
        
        # Get the time value 
        T = T + dt
        T_array = np.append(T_array, T)
        
        for j in range(1,n):
            P_Euler_new[j,0] = 0.5*(P_Euler[j-1,0]+P_Euler[j+1,0]) - dt/(2*dx)*(Q_Euler[j+1,0]-Q_Euler[j-1,0])
            Q_Euler_new[j,0] = 0.5*(Q_Euler[j-1,0]+Q_Euler[j+1,0]) - dt/(2*dx)*a_square*(P_Euler[j+1,0]-P_Euler[j-1,0]) - dt*teta/4*(Q_Euler[j-1,0]*abs(Q_Euler[j-1,0])/P_Euler[j-1,0] + Q_Euler[j+1,0]*abs(Q_Euler[j+1,0])/P_Euler[j+1,0]) 
            #0.5*(Q[j-1]+Q[j+1]) - dt/(2*dx)*a_square*(P[j+1,t]-P[j-1,t]) - dt*Lambda/(4*D)*(Q[j-1,t]*abs(Q[j-1,t])/P[j-1,t] + Q[j+1,t]*abs(Q[j+1,t])/P[j+1,t])
        
        #Kompatibilitätsbedingungen
        P_Euler_new[0,0] = 4.75*T*T+1 #2.7067*exp(T) #node 0
        P_Euler_new[n,0] = P_Euler[n,0]         #last node, Kompatibilität

        Q_Euler_new[0,0] = 10*T*T #5.4134*exp(T) #node 0
        Q_Euler_new[n,0] = Q_Euler[n,0]          #last node

        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        P_Matrix = np.append(P_Matrix, P_Euler, axis = 1)
        Q_Matrix = np.append(Q_Matrix, Q_Euler, axis = 1)
        
        # print every 300 general steps
        GeneralSteps = GeneralSteps + 1
        
        if GeneralSteps%50 == 0:
            print("General Step:" + str(GeneralSteps) + " ,Time: " + str(T))
            
    return dt, GeneralSteps, X_array[:-1], T_array, P_Matrix, Q_Matrix


# Euler Gleichung mit Simple Upwind
def EulerGleichung_SimpleUp():
    # initial conditions
    P_Euler = np.zeros((n+1,1)) #Spaltenvektor
    P_Euler_new = np.zeros((n+1,1)) #Spaltenvektor
    
    x = 0
    X_array = [x]
    for i in range(0,n+1):
        P_Euler[i,0] = exp(beta*x) #Spaltenvektor
        x = x + dx
        X_array = np.append(X_array,x)
        
    Q_Euler = np.zeros((n+1,1))  #initialisiere
    Q_Euler_new = np.zeros((n+1,1))      
            
    P_Matrix = P_Euler
    Q_Matrix = Q_Euler
    
    T = 0
    GeneralSteps = 0
    T_array = [T]
    
    while T < 2:
        
        # Q_Euler/P_Euler computes the velocity
        # Calculate dt with the CFL condition
        #dt = (CFL*dx)/(np.amax(Q_Euler/P_Euler + a_square*np.ones((1,n+1))))
        dt = 0.0001
        
        # Get the time value 
        T = T + dt
        T_array = np.append(T_array, T)
        
        for j in range(1,n+1):
            P_Euler_new[j,0] = P_Euler[j,0] + dt/dx*(Q_Euler[j-1,0]-Q_Euler[j,0])
            Q_Euler_new[j,0] = Q_Euler[j,0] + dt/dx*((Q_Euler[j-1,0]*Q_Euler[j-1,0])/(P_Euler[j-1,0]) + a_square*P_Euler[j-1,0] - (Q_Euler[j,0]*Q_Euler[j,0])/(P_Euler[j,0])-a_square*P_Euler[j,0]) - dt*teta/2*(Q_Euler[j,0]*abs(Q_Euler[j,0])/P_Euler[j,0]) 
            #0.5*(Q[j-1]+Q[j+1]) - dt/(2*dx)*a_square*(P[j+1,t]-P[j-1,t]) - dt*Lambda/(4*D)*(Q[j-1,t]*abs(Q[j-1,t])/P[j-1,t] + Q[j+1,t]*abs(Q[j+1,t])/P[j+1,t])
        
        #Kompatibilitätsbedingungen
        P_Euler_new[0,0] = 4.75*T*T+1 #2.7067*exp(T) #node 0
        #P_Euler_new[n,0] = P_Euler[n,0]         #last node, Kompatibilität

        Q_Euler_new[0,0] = 10*T*T #5.4134*exp(T) #node 0
        #Q_Euler_new[n,0] = Q_Euler[n,0]          #last node

        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        P_Matrix = np.append(P_Matrix, P_Euler, axis = 1)
        Q_Matrix = np.append(Q_Matrix, Q_Euler, axis = 1)
        
        # print every 300 general steps
        GeneralSteps = GeneralSteps + 1
        
        if GeneralSteps%50 == 0:
            print("General Step:" + str(GeneralSteps) + " ,Time: " + str(T))

    return dt, GeneralSteps, X_array[:-1], T_array, P_Matrix, Q_Matrix


# Euler Gleichung mit Simple Upwind
def Weymouth_SimpleUp():
    # initial conditions
    P_Euler = np.zeros((n+1,1)) #Spaltenvektor
    P_Euler_new = np.zeros((n+1,1)) #Spaltenvektor
    
    x = 0
    X_array = [x]
    for i in range(0,n+1):
        P_Euler[i,0] = exp(beta*x) #Spaltenvektor
        x = x + dx
        X_array = np.append(X_array,x)
        
    Q_Euler = np.zeros((n+1,1))  #initialisiere
    Q_Euler_new = np.zeros((n+1,1))      
            
    P_Matrix = P_Euler
    Q_Matrix = Q_Euler
    
    T = 0
    GeneralSteps = 0
    T_array = [T]
    
    # in order to fulfill CFL condition, we need 
    # dt*a <= dx -> dt*1 <= 0.0625 -> dt = 0.05
    dt = 0.001
    
    while T < 2:
        
        # Get the time value 
        T = T + dt
        T_array = np.append(T_array, T)
        
        for j in range(1,n+1):
            P_Euler_new[j,0] = P_Euler[j,0] - dt/dx*(Q_Euler[j,0]-Q_Euler[j-1,0])
            Q_Euler_new[j,0] = Q_Euler[j,0] - dt/dx*a_square*(P_Euler[j,0]-P_Euler[j-1,0]) - dt*teta/2*(Q_Euler[j,0]*abs(Q_Euler[j,0])/P_Euler[j,0]) 
            #0.5*(Q[j-1]+Q[j+1]) - dt/(2*dx)*a_square*(P[j+1,t]-P[j-1,t]) - dt*Lambda/(4*D)*(Q[j-1,t]*abs(Q[j-1,t])/P[j-1,t] + Q[j+1,t]*abs(Q[j+1,t])/P[j+1,t])
        
        #Kompatibilitätsbedingungen
        P_Euler_new[0,0] = 4.75*T*T+1 #2.7067*exp(T) #node 0
        #P_Euler_new[n,0] = P_Euler[n,0]         #last node, Kompatibilität

        Q_Euler_new[0,0] = 10*T*T #5.4134*exp(T) #node 0
        #Q_Euler_new[n,0] = Q_Euler[n,0]          #last node

        P_Euler = P_Euler_new
        Q_Euler = Q_Euler_new
        
        P_Matrix = np.append(P_Matrix, P_Euler, axis = 1)
        Q_Matrix = np.append(Q_Matrix, Q_Euler, axis = 1)
        
        # print every 300 general steps
        GeneralSteps = GeneralSteps + 1
        
        if GeneralSteps%50 == 0:
            print("General Step:" + str(GeneralSteps) + " ,Time: " + str(T))
            
    return dt, GeneralSteps, X_array[:-1], T_array, P_Matrix, Q_Matrix


# Euler Gleichung mit LxF
if __name__ == '__main__':
    
    EulerPlt = 1
    WeymouthBlt = 0
    SimpleUp_Euler_Blt = 0
    SimpleUp_Weymouth_Blt = 0
    
    if (EulerPlt == 1):
        # Plot Euler
        dt1, steps1, X_Arr, T_Arr, P_Matrix, Q_Matrix = EulerGleichung_LxF()
        X_mesh, T_mesh = np.meshgrid(X_Arr, T_Arr)
        P_Matrix_T = P_Matrix.transpose()
        Q_Matrix_T = Q_Matrix.transpose()
        
        plt.figure(1)
        ax = plt.axes(projection='3d')
        ax.dist = 11
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_zlabel('p(t,x)')
        ax.locator_params(axis='x', nbins = 6)
        ax.locator_params(axis='y', nbins = 5)
        ax.locator_params(axis='z', nbins = 5)
        ax.plot_surface(X_mesh, T_mesh, P_Matrix_T)#, rstride=1, cstride=1) ###!
                    #cmap='viridis', edgecolor='none')
        title = 'p(t,x) with Euler Equation with step size: ' + str(dt1)
        ax.set_title(title)
        title_file1 = 'Euler_P_' + str(dt1) + '.png'
        plt.savefig(title_file1, dpi=300)
        
        plt.figure(2)
        ax2 = plt.axes(projection='3d')
        ax2.dist = 11
        ax2.set_xlabel('x')
        ax2.set_ylabel('t')
        ax2.set_zlabel('q(t,x)')
        ax2.locator_params(axis='x', nbins = 6)
        ax2.locator_params(axis='y', nbins = 5)
        ax2.locator_params(axis='z', nbins = 5)
        ax2.plot_surface(X_mesh, T_mesh, Q_Matrix_T) #, rstride=1, cstride=1,
                    #cmap='viridis', edgecolor='none')
        title = 'q(t,x) with Euler Equation with step size: ' + str(dt1)
        ax2.set_title(title)
        title_file2 = 'Euler_Q_' + str(dt1) + '.png'
        plt.savefig(title_file2, dpi=300)

    if (WeymouthBlt == 1):
    # Plot Weymouth
        dt2, steps2, X_Arr2, T_Arr2, P_Matrix2, Q_Matrix2 = Weymouth_LxF()
        X_mesh2, T_mesh2 = np.meshgrid(X_Arr2, T_Arr2)
        P_Matrix_T2 = P_Matrix2.transpose()
        Q_Matrix_T2 = Q_Matrix2.transpose()
        
        plt.figure(3)
        ax3 = plt.axes(projection='3d')
        ax3.dist = 11
        ax3.set_xlabel('x')
        ax3.set_ylabel('t')
        ax3.set_zlabel('p(t,x)')
        ax3.locator_params(axis='x', nbins = 6)
        ax3.locator_params(axis='y', nbins = 5)
        ax3.locator_params(axis='z', nbins = 5)
        ax3.plot_surface(X_mesh2, T_mesh2, P_Matrix_T2)#, rstride=1, cstride=1) ###!
                    #cmap='viridis', edgecolor='none')
        title = 'p(t,x) with Weymouth Equation with step size: ' + str(dt2)
        ax3.set_title(title)
        title_file3 = 'Weymouth_P_' + str(dt2) + '.png'
        plt.savefig(title_file3, dpi=300)
        
        plt.figure(4)
        ax4 = plt.axes(projection='3d')
        ax4.dist = 11
        ax4.set_xlabel('x')
        ax4.set_ylabel('t')
        ax4.locator_params(axis='x', nbins = 6)
        ax4.locator_params(axis='y', nbins = 5)
        ax4.locator_params(axis='z', nbins = 5)
        ax4.set_zlabel('q(t,x)')
        ax4.plot_surface(X_mesh2, T_mesh2, Q_Matrix_T2) #, rstride=1, cstride=1,
                    #map='viridis', edgecolor='none')
        title = 'q(t,x) with Weymouth Equation with step size: ' + str(dt2)
        ax4.set_title(title)
        title_file4 = 'Weymouth_Q_' + str(dt2) + '.png'
        plt.savefig(title_file4, dpi=300)
        
    if (SimpleUp_Euler_Blt == 1):
        dt3, steps3, X_Arr3, T_Arr3, P_Matrix3, Q_Matrix3 = EulerGleichung_SimpleUp()
        X_mesh3, T_mesh3 = np.meshgrid(X_Arr3, T_Arr3)
        P_Matrix_T3 = P_Matrix3.transpose()
        Q_Matrix_T3 = Q_Matrix3.transpose()
        
        plt.figure(5)
        ax5 = plt.axes(projection='3d')
        ax5.dist = 11
        ax5.set_xlabel('x')
        ax5.set_ylabel('t')
        ax5.set_zlabel('p(t,x)')
        ax5.locator_params(axis='x', nbins = 6)
        ax5.locator_params(axis='y', nbins = 5)
        ax5.locator_params(axis='z', nbins = 5)
        ax5.plot_surface(X_mesh3, T_mesh3, P_Matrix_T3)#, rstride=1, cstride=1) ###!
                    #cmap='viridis', edgecolor='none')
        title = 'p(t,x) Simple Upwind Euler Equation, step size: ' + str(dt3)
        ax5.set_title(title)
        title_file5 = 'SimpleUp_Euler_P_' + str(dt3) + '.png'
        plt.savefig(title_file5, dpi=300)
        
        plt.figure(6)
        ax6 = plt.axes(projection='3d')
        ax6.dist = 11
        ax6.set_xlabel('x')
        ax6.set_ylabel('t')
        ax6.locator_params(axis='x', nbins = 6)
        ax6.locator_params(axis='y', nbins = 5)
        ax6.locator_params(axis='z', nbins = 5)
        ax6.set_zlabel('q(t,x)')
        ax6.plot_surface(X_mesh3, T_mesh3, Q_Matrix_T3) #, rstride=1, cstride=1,
                    #map='viridis', edgecolor='none')
        title = 'q(t,x) Simple Upwind Euler Equation, step size:  ' + str(dt3)
        ax6.set_title(title)
        title_file6 = 'SimpleUp_Euler_Q' + str(dt3) + '.png'
        plt.savefig(title_file6, dpi=300)
        
        
    if (SimpleUp_Weymouth_Blt == 1):
        dt4, steps4, X_Arr4, T_Arr4, P_Matrix4, Q_Matrix4 = Weymouth_SimpleUp()
        X_mesh4, T_mesh4 = np.meshgrid(X_Arr4, T_Arr4)
        P_Matrix_T4 = P_Matrix4.transpose()
        Q_Matrix_T4 = Q_Matrix4.transpose()
        
        plt.figure(7)
        ax7 = plt.axes(projection='3d')
        ax7.dist = 11
        ax7.set_xlabel('x')
        ax7.set_ylabel('t')
        ax7.set_zlabel('p(t,x)')
        ax7.locator_params(axis='x', nbins = 6)
        ax7.locator_params(axis='y', nbins = 5)
        ax7.locator_params(axis='z', nbins = 5)
        ax7.plot_surface(X_mesh4, T_mesh4, P_Matrix_T4)#, rstride=1, cstride=1) ###!
                    #cmap='viridis', edgecolor='none')
        title = 'p(t,x) Simple Upwind Weymouth Equation, step size: ' + str(dt4)
        ax7.set_title(title)
        title_file7 = 'SimpleUp_Euler_P_' + str(dt4) + '.png'
        plt.savefig(title_file7, dpi=300)
        
        plt.figure(8)
        ax8 = plt.axes(projection='3d')
        ax8.dist = 11
        ax8.set_xlabel('x')
        ax8.set_ylabel('t')
        ax8.locator_params(axis='x', nbins = 6)
        ax8.locator_params(axis='y', nbins = 5)
        ax8.locator_params(axis='z', nbins = 5)
        ax8.set_zlabel('q(t,x)')
        ax8.plot_surface(X_mesh4, T_mesh4, Q_Matrix_T4) #, rstride=1, cstride=1,
                    #map='viridis', edgecolor='none')
        title = 'q(t,x) Simple Upwind Weymouth Equation, step size:  ' + str(dt4)
        ax8.set_title(title)
        title_file8 = 'SimpleUp_Euler_Q' + str(dt4) + '.png'
        plt.savefig(title_file8, dpi=300)