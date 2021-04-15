#Aidan Walk
#ASTR 260-001
#13 April 2021, 17.00
#Homework 10 - More Fun With ODE's

#Problem 1
import numpy as np

def secondDerivative(x, t):
    ''' d^2x/dt^2 - (dx/dt)^2 + x + 5 = 0 
    '''
    #let y = dx/dt -> dy/dt = d^2x/dt^2
    y = firstDerivative(x, t)
    dyDt = y^2 - x - 5
    
    return dyDt

def firstDerivative(x, t):
    '''
    '''
    
    dxDt = np.sqrt(secondDerivVal + x + 5)

    return dxDt
    
def leapFrog(x0, dxDt0, times=None, deriv=None):
    '''Computes the solution of a second order ODE via the leapfrog method
    '''
    #set current xVal to initial xVal
    #set current vVal to initial vVal
    for t in times:
        vAt_tPlusHalfH = v(t) + 1/2*h*f(x(t), t)
        xAt_tPlusH = x(t) + h*vAt_tPlusHalfH
        
        k = h*f(xAt_tPlusH, t+h)
        
        vAt_tPlusH = vAt_tPlusHalfH + 1/2*k
        vAt_tPlus3HalfH = vAt_tPlusHalfH + k
       
    pass
    
def RK4(x0, t0, h=None, deriv=None):
    '''Computes the first half value for the leap frog method using RK4
       x0 = initial value
       h = step size
       deriv = derivative we are solving'''
    t = t0
    
    k1 = h * deriv(x0, t)
    k2 = h * deriv(x0 + 1/2*k1, t + h/2)
    k3 = h * deriv(x0 + 1/2*k2, t + h/2)
    k4 = h * deriv(x0 + k3, t + h)
    
    nextValue = x0 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    
    return nextValue
    
if __name__ == '__main__':
    #time intervals
    tNaught = 0
    tFinal = 50 
    step = 0.001
    times = np.arange(tNaught, tFinal, step)
    
    #initial conditions
    xNaught = 1
    dxDt = 0
    
    
    #compute first half value using RK4
    x_hHalf_2ndDeriv = RK4(xNaught, tNaught,
                           h=step, deriv=secondDerivative)
    x_hHalf = RK4(x_hHalf_2ndDeriv, tNaught,
                  h=step, deriv=firstDerivative)
    
    #solve second derivative and pass into first
    
    
    #solve first derivative to obtain final value
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    