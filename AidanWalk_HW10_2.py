#Aidan Walk
#ASTR 260-001
#13 April 2021, 17.00
#Homework 10 - More Fun With ODE's

#Problem 2
import numpy as np
import matplotlib.pyplot as plt
import time


#Fix this
G = 6.6743e-11 #kg^-1 m^3 s^-2 - Gravitational Constant
M = 1.9885e30 #kg - Mass of the sun
AU = 1.496e11 #1 astronomical unit in meters

def f(r, t=np.nan):
    """derivative function to pass to rk4
    pass state vector r = [x, y, vx, vy]"""
    x, y, vx, vy = r
    rcubed = np.sqrt(x**2 + y**2)**3

    fx = vx #return 1st parameter
    fy = vy #return 2nd parameter
    fvx = -G*M*x/rcubed #return 3rd parameter
    fvy = -G*M*y/rcubed #return last parameter
    
    return np.array([fx, fy, fvx, fvy])

def rk4_step(r=None, h=None, f=None, t=np.nan):
    '''computes k1 through k4 for RK4 method
       returns next value (1/6 * (k1 + 2*k2 + 2*k3 + k4))'''
    
    k1 = h*f(r, t)
    k2 = h*f(r + 0.5*k1, t + 0.5*h)
    k3 = h*f(r + 0.5*k2, t + 0.5*h)
    k4 = h*f(r + k3, t + h)
    
    return 1.0/6*(k1 + 2*k2 + 2*k3 + k4)

def run_rk4_fixed(initial_state=None, initial_h=None, tmax=None):
    '''computes RK4 with a fixed step size'''
    t=0
    times = np.arange(t, tmax, initial_h) #times to loop over
    
    stateVector = stateVector_array = initial_state #initialize array and state vec
    for t in times:
        stateVector = stateVector + rk4_step(r=stateVector,
                               h=initial_h,
                               f=f)
        stateVector_array = np.vstack((stateVector_array, stateVector)) #append data 
    
    xpos, ypos = stateVector_array[:,0], stateVector_array[:,1]
    return xpos, ypos
    
def run_rk4_adaptive(initial_state=None, initial_h=None, tmax=None):
    r = initial_state
    h = initial_h
    xpoints = []
    ypoints = []
    
    t=0
    while t<tmax:
        #Do one large step
        r1 = r + rk4_step(r=r, h=2*h, f=f)
        #Do two small steps
        r2 = r + rk4_step(r=r, h=h, f=f)
        r2 = r2 + rk4_step(r=r2,h=h, f=f)

        #calculate value of rho
        ex = 1.0/30*(r[0]-r1[0])
        ey = 1.0/30*(r[1]-r1[1])
        rho = 30*h*delta/np.sqrt(ex**2+ey**2)

        #calculate new values of t, h, r
        #update points if appropriate
        if rho>=1.0:
            h = 2*h
            print('h =', h)
            t += h
            r = r2
            xpoints.append(r[0])
            ypoints.append(r[1])
            
        else:
            h = h*rho**(1/4)
            r = r
            print(h)
            
    return np.array([xpoints, ypoints])

def split_stateVec(x):
    '''
    '''
    pos, vel = np.split(x, 2)
    return pos, vel

if __name__ == "__main__":
    h0_adaptive = 1.0e5 #initial step size for RK4 adaptive
    h_fixed = 3.0e4 #step size to accurately model orbit with fixed step size
    tmax_adaptive = 1.0e9 #total time for adaptive Rk4
    tmax_fixed = 3e9 #total time for fixed RK4

    delta = 1e3/(365.25*24*3600) #meters accuracy per second

    x0, y0 = 4e12, 0 #starting pos, 4 billion kilometers
    vx0, vy0 = 0, 500 #starting velocity, m/s
    r0 = np.array([x0, y0, vx0, vy0]) #<-standard layout for all statevectors in this code

    #calculate RK4 fixed
    print('Calculating two cometary orbits using fixed step size...(I really need to optimize this)\n')
    print('For decent accuracy, h =', h_fixed, 'meters') 
    #time calculation 
    start_time = time.time()
    xpos_RK4Fixed, ypos_RK4Fixed = run_rk4_fixed(initial_state=r0,
                                                 initial_h=h_fixed, 
                                                 tmax=tmax_fixed)
                                                 
    end_time = time.time() #stop time
    print('--- Calculation duration: %s seconds ---' % (end_time - start_time))
    print('As h gets bigger, the energy lost in each orbit becomes greater.')
    print('Any larger of a step size results in the orbit spiraling into the sun')
    
    #calculate adaptive RK4
    #CANNOT GET THIS TO WORK
    #RK4_adaptive = run_rk4_adaptive(initial_state=r0, initial_h=h0_adaptive, tmax=tmax_adaptive)
    
    
    #Make the plot
    #Plot RK4 fixed
    plt.plot(xpos_RK4Fixed/AU, ypos_RK4Fixed/AU, alpha = 0.5)
    
    #plot RK4 adaptive
    #plt.plot(stateVec_RK4Fixed[:,0]/AU, stateVec_RK4Fixed[:,1]/AU, alpha = 0.5)
    #plt.plot(stateVec_RK4Fixed[:,0]/AU, stateVec_RK4Fixed[:,1]/AU, 'k.')
    plt.show()
