#By Moises Sanchez (inspired by Burden's Numerical Analysis)

#This is an Euler's Method First Order ODE Solver
#let a,b be our endpoints, N be an integer, and o be our intial condtion
#if you want to find actual error, you can put the general solution to the ODE at true_f(t)
#if you are not looking for the actual error then just comment out lines 13,14,25,26,28,29

import numpy as np

#place ODE that you want to solve
def f(w,t):
    return (w - t ** (2) + 1)
#given the actual solution to the ODE, plug here
def true_f(t):
    return (t+1) ** 2 - 0.5 * np.exp(1) ** t
    
#Euler's method algorithum
def Euler_Method_Ode(a,b,N,o):
    h = (b-a) / N
    t = a
    w = o
    print('running process ...')
    for k in range(1,N+1):
        w = w + h * f(w,t)
        t = a + (k * h)
        y_i = true_f(t)
        actual_error = np.abs(y_i-w)
    print('Our approximate value of our ODE is ' + str(w) + ' at t = ' + str(t))
    print('Our actual value of the ODE is ' + str(y_i) + ' at t = ' + str(t))
    print('Our actual error is ' + str(actual_error))
    
