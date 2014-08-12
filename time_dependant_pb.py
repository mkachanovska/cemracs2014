#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import os
from time import time
import numpy as np
from pylab import *
import scipy.io
from scipy.special import airy
import copy

#---------------CONSTANTS--------------------#
# mesh size
N      = 1024
# domain ]-L ; H[, (mL=-L)
mL     = -1
H      =  10
# space step 
dx     = (H-mL)/N
print "dx", dx
# time step
dt     = 0.5*dx
print "dt",dt
# Number of time steps
Ntime = 50
# constants 
me     = 1#9.11e-31               #electron mass (kg)
e      = 1#1.6e-19               #electron charge (coulombs)
eps0   = 1#8.85e-12               #permeability  (F/m)
nu     = 0                     #friction
B0     = 0                      #given
c      = 1#3e8 



#alpha(-L)
alphamL = 2

#-------------airy functions------------------#
#return A
#def airy_one(x):
#    [A,B,C,D] = airy(x)
#    return A
#return A'
#def airy_two(x):
#    [A,B,C,D] = airy(x)
#    return B
#-----------------------------------------------#

#-------------------------- NE --------------------------#
# Ne constants definition
wc    = e*abs(B0)/me
print "omega_c" ,wc
# Ne function
def Ne(x) :  
    a     = wc/(c*(wc*wc-c*c))
    b     = 1/(wc*wc-c*c)
    kappa = e*e/(eps0*me)

    a2    = (a*a-b*b)*kappa*kappa
    a1    = (2*kappa*b + b*x*kappa) 
    a0    = 1 + x
    
    #D     = a1*a1 + 4 * a2 * a0
    
    return  - a1 /(2*a2)

#NE    = map(lambda x: 3*me/(e*e),  X)
omega = np.sqrt(7)
wp = np.sqrt(3)
gamma = np.sqrt(omega*omega-wp*wp)
#---------------------------------------------------------#

#--------------G for left BC -------------------#
alpha = 2

def g(t) : 
    return np.exp(complex(0,omega*t))*np.exp(-gamma*mL)*(-gamma+complex(0,alphamL))
#-----------------------------------------------#
#------table and functions initialisation-------#
#-----------------------------------------------#
X    = map(lambda i: mL + i*dx, range(N+1))
NE   = map(lambda x: 3*me/(e*e),  X)
X12  = map(lambda i: mL+0.5*dx + i*dx, range(N))
NEy  = map(lambda x: 3*me/(e*e),  X12)
Ex   = map(lambda x: 0, X)
ux   = map(lambda x: 0, X)
H    = map(lambda x: (gamma/complex(0,omega))*np.exp(-gamma*x)*np.exp(complex(0,-omega*dt/2)), X)
Ey   = map(lambda x: np.exp(-gamma*x),X12)
uy   = map(lambda x: -e/(me*complex(0,omega))*np.exp(-gamma*x), X12)

#H12 = map(lambda x: -np.exp(-alpha*dt/2)*np.exp(-alpha*x),X)

tux = copy.deepcopy(ux)
tuy = copy.deepcopy(uy)
tEx = copy.deepcopy(Ex)
tEy = copy.deepcopy(Ey)
H12 = copy.deepcopy(H)
scipy.io.savemat('H120.mat', {'H120':H12})



#---------------------------------------------------------#

#--------------------K1&K2 s def function-------------------#
def Kcoeff(x): 
    K1     =  1 + nu * dt/2 + dt*dt*x*e*e/(4*me)
    K2     =  1 - nu * dt/2 - dt*dt*x*e*e/(4*me)
    #K1x    =  K1 + beta*beta/(4*K1)
    #K2x    =  K2 - beta*beta/(4*K1)
    K1x = 0
    K2x = 0
    return [K1,K2,K1x,K2x]
#-----------------------------------------------------------#
#------------------------time loop--------------------------#
#-----------------------------------------------------------#


t = 1

#coeff alpha (1/alpha_Ey)Ey in left BC
alphaey = complex(0,alphamL)/2 - 1/dx
print 'alphaey', alphaey

while (t<=Ntime):
    print "iter : ",t,"time : ", t*dt

    #---------------- H -> H12
    H12[0] = H[0] -  dt/dx * (Ey[0] \
                                    - (1/alphaey) * (g((t-1)*dt) -  complex(0,alphamL)/2 * Ey[0] \
                                                         - (1/dx)* Ey[0]))

    i = 1
    while (i<N):
        H12[i] = H[i] - dt/dx*(Ey[i] - Ey[i-1])
        i      = i +1
    H12[N] = H[N]

    #------------------- ux -> tux
    #for i in range(N-1):
        
    #   # [K1,K2,K1x,K2x] = Kcoeff(NE[i])
    #    tux[i] =  0 #(1/K1x) * (K2x*ux[i] + gamma*Ex[i] + ((beta*gamma)/(2*K1))*Ey[i]\
    #                #           - beta*gamma*dt/(4*K1)*((H12[i] - H12[i-1])/dx)   \
    #                #           + (beta/2)*(K2/K1 + 1)*uy[i])
        
    #left BC (if i = N) (H(i+1) = 0) ????
    #[K1,K2,K1x,K2x] = Kcoeff(NE[N+1])
    #tux[N] =  0#(1/K1x)*(K2x*ux[N+1] + gamma*Ex[N+1] + ((beta*gamma)/(2*K1))*Ey[N]\
                 #             - beta*gamma*dt/(4*K1)*(H12[N+1] - H12[N])/dx   \
                 #             + (beta/2)*(K2/K1 + 1)*uy[N])

    scipy.io.savemat('uy.mat',  {'uy': uy}); 

    #-------------------- u y->tuy

    for i in range(N):
        [K1,K2,K1x,K2x] = Kcoeff(NE[i])
        tuy[i] = (1/K1) * (K2 * uy[i] -dt * e / me * Ey[i] +dt*dt*e/(2*me) * (H12[i+1] - H12[i])/dx)

    #------------------- E -> tE
                        

    for i in range(N):
       # tEx[i] = Ex[i] + e*NE[i]*dt* (tux[i] + ux[i])/2

        tEy[i] = Ey[i] - (dt) * (H12[i+1] - H12[i])/dx   #+(dt*e*NEy[i]/(2*me))*(tuy[i] + uy[i])

    #left BC (if i = N)(H(i+1) = 0) ????
  
    #tEx[N-1]    = Ex[N-1] - (dt*e*NE[N-1])*(tux[N-1] + ux[N-1])/2
    #tEx[N]  = Ex[N] - (dt*e*NE[N])*(tux[N] + ux[N])/2
    #scipy.io.savemat('Ey.mat',  {'arr': tEy});
     
    

    #----------f^n <- f^n+1
    ux = copy.deepcopy(ux)
    Ex = copy.deepcopy(tEx)
    uy = copy.deepcopy(tuy)
    Ey = copy.deepcopy(tEy)
    H  = copy.deepcopy(H12)

    time = t*dt

    scipy.io.savemat('H12.mat',  {'H12': H12})
    scipy.io.savemat('Eyo.mat',  {'Eyo': Ey})
    scipy.io.savemat('tuy.mat',  {'tuy': tuy})
    Eyex   = map(lambda x: np.exp(complex(0,omega*time))*np.exp(-gamma*x),X12)
    Hex    = map(lambda x: (gamma/complex(0,omega))*np.exp(-gamma*x)*np.exp(complex(0,omega*(time-dt/2))) ,X)
    Uyex   = map(lambda x: -e/(me*complex(0,omega))*np.exp(complex(0,omega*time)-gamma*x),X12)

    scipy.io.savemat('Hex.mat',  {'Hex': Hex})
    scipy.io.savemat('Eyex.mat',  {'Eyex': Eyex})
    scipy.io.savemat('uyex.mat',  {'uyex': Uyex})   

    temp = 0
    for i in range(N):
        temp = temp + np.abs((H[i]-Hex[i]))**2
    print time, 'norm err H', np.sqrt(dx*temp)

    temp = 0
    for i in range(N) : 
        temp = temp + np.abs((uy[i]-Uyex[i]))**2
    print 'norm err uy',np.sqrt(dx*temp)
    
    temp = 0
    for i in range(N) : 
        temp = temp + np.abs((Ey[i]-Eyex[i]))**2
    print 'norm err Ey', np.sqrt(dx*temp)                         

    t = t+1           



