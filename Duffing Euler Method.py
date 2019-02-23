#By Moises Jornales Sanchez and Eric Barrios
#Subject: Math 104B
#This is an ODE solver using the Euler's Method for the Duffing Oscillator
#To change the values for the ODE check the values in run_Euler, and runner_interval
#---------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import astropy.table as ast
import time


#this sets up to what value we are aiming for given our steps
def interval(a,h,N):
    for i in range(N+1):
        c_int = a + i*h
    return(c_int)

#the actual Euler's Method
def Euler_Method(f_1,f_2,energy_equation,a,h,N,IV_1,IV_2,IPos,IVel,delta,alpha,beta,omega,gamma):
    start_time = time.time()

    steps = runner_interval()
    t = np.arange(a,(steps+h)/(2*np.pi/omega),h) #creates our mesh
    t = np.round(t,decimals = 3) #stops mesh from giving unwanted recurring decimals
    w_1 = np.zeros(len(t)) #initalizes our w values for approx.
    w_2 = np.zeros(len(t)) #initalizes our w values for approx.
    
    #setting our intial conditions
    t[0],w_2[0] = IV_1 #please change w_2 for diff starting
    t[0],w_1[0] = IV_2 #please change w_2 for diff starting
    variance_energy = 0
    mean_energy = 0
    
    for k in range(1,len(t)):
        w_2[k] = w_2[k-1] + (h)*f_2(t[k-1],w_2[k-1],w_1[k-1],delta,alpha,beta,omega) 
        w_1[k] = w_1[k-1] + (h)*w_2[k]
        
    if ((gamma == 0) and (delta == 0)): #given we have a non-driven undampened oscillator we calc. mean & variance energy
       
        for k in range(1,len(t)): #this calculates the mean energy
            mean_energy += (energy_equation(w_1,w_2,beta,alpha))
        mean_energy = sum(mean_energy) / len(t)
        
        for k in  range(1,len(t)):
            variance_energy += ((energy_equation(w_1,w_2,beta,alpha)) - mean_energy) ** 2
        variance_energy = sum(variance_energy) / len(t)
        
        
    print("--- %s seconds ---" % (time.time() - start_time))
    print("Average of the energy for this system is: " + str(mean_energy))
    print("Variance of the energy for this system is: " + str(variance_energy)) 
    plt.figure(1)
    plt.plot(w_1,w_2, color = "b", label = 'Trajectory')
    plt.title("Duffing Phase Diagram: Euler's Method, N = " +str(len(t)-1))
    plt.xlabel('postion')
    plt.ylabel('velocity')
    plt.legend(loc =4)
    plt.grid()
    plt.savefig('C:/Users/moijs/Desktop/Math 104B/Euler_TimeSeries.png',fmt='PNG',dpi =100)
    plt.figure(2)
    plt.plot(t,w_1, color = "r")
    plt.title("Duffing Time Trace: Euler's Method, N = " +str(len(t)-1))
    plt.xlabel('frequency')
    plt.ylabel('postion')
    plt.grid()
    plt.savefig('C:/Users/moijs/Desktop/Math 104B/EulerMod_Phase.png',fmt='PNG',dpi =100)
    return (ast.Table([t,w_1,w_2],names = ['t_i','Euler: Velocity','Euler: Postion']))
      
#these run the ODE Solver
#------------------------------------------------

def run_Euler():
    f_1 = lambda t,w_1: w_1 #the postion of the velocity of the ODE
    f_2 = lambda t,w_1,w_2,delta,alpha,beta,omega: -delta * w_1 - alpha * w_2 - beta * w_2 ** 3 + gamma * np.cos(omega * t)  #the acc. of the ODE
    energy_equation = lambda w_1,w_2,beta,alpha:  w_1  ** 2 /2 - alpha * w_2 ** 2 /2 - beta * w_2 ** 4 / 4 
    delta = 0.3 #damping constant
    alpha = -1 #linear stiffness
    beta = 1 #second damping constant
    gamma = 0.5 #amplitude
    omega = 1.2 #frequency
    a = 0.0 #initial time
    h = 0.025 #steps
    N = 100000 #interations
    IPos = 1 #intital postion
    IVel = 0 #initial velocity
    IV_1 = (a,IPos) #initial postion and time
    IV_2 = (a,IVel) #initial velocity and time
    return (Euler_Method(f_1,f_2,energy_equation,a,h,N,IV_1,IV_2,IPos,IVel,delta,alpha,beta,omega,gamma))

def runner_interval():
    a= 0.0 #initial time
    h = 0.025 #steps
    N = 100000 #initial postion
    return(interval(a,h,N))

run_Euler()