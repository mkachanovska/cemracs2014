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
N      = 1000
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
Ntime = 1
# constants 
me     = 1#9.11e-31               #electron mass (kg)
e      = 1#1.6e-19               #electron charge (coulombs)
eps0   = 1#8.85e-12               #permeability  (F/m)
nu     = 0                    #friction
B0     = 0                      #given
c      = 1#3e8 



#alpha(-L)
alphamL = 2

# Ne constants definition
wc    = e*abs(B0)/me
print "omega_c" ,wc

wp = np.sqrt(3);
omega=wp;
gamma = 2
#--------------G for left BC -------------------#
def g(t) :
    return (-gamma + complex(0,alphamL))*np.exp(complex(-gamma*mL,omega*t))

#-----------------------------------------------#
#------table and functions initialisation-------#
#-----------------------------------------------#
X    = map(lambda i: mL + i*dx, range(N+1))
NE   = map(lambda x: wp**2*me/(e**2),  X)
X12  = map(lambda i: mL+0.5*dx + i*dx, range(N))
NEy  = map(lambda x: wp**2*me/(e**2),  X12)
Ex   = map(lambda x: np.exp(-gamma*x), X)
ux   = map(lambda x: -(e/me)*(1/(complex(0,omega)))*x, Ex)
H    = map(lambda x: 0,X)#(gamma/complex(0,omega))*np.exp(-gamma*x), X)
Ey   = map(lambda x: 0,X12)#np.exp(-gamma*x),X12)
uy   = map(lambda x: 0,X12)#-(e/me)*(1/complex(0,omega))*np.exp(-gamma*x), X12)

tux = copy.deepcopy(ux)
tuy = copy.deepcopy(uy)
tEx = copy.deepcopy(Ex)
tEy = copy.deepcopy(Ey)
H12 = copy.deepcopy(H)
scipy.io.savemat('H120.mat', {'H120':H12})

#-----------------------------------------------------------#
#------------------------time loop--------------------------#
#-----------------------------------------------------------#
def Kcoeff(x): 
    #K1     =  1 + nu * dt/2 + dt*dt*x*e*e/(4*me)
    #K2     =  1 - nu * dt/2 - dt*dt*x*e*e/(4*me)
    K1x =  1 + nu*dt/2+(dt**2)*e*x/(4*me);
    K2x =  1 - nu*dt/2-(dt**2)*e*x/(4*me);
    K1  = 0
    K2  = 0
    return [K1,K2,K1x,K2x]
[K1,K2,K1x,K2x] = Kcoeff(wp**2*me/e**2)


#coeff alpha (1/alpha_Ey)Ey in left BC
alphaey = complex(0,alphamL)/2 - 1/dx
print 'alphaey', alphaey
t = 1
Ntime = 120
while (t<=Ntime):
    time = t*dt

    
    for i in range(N+1):
        tux[i] =  (1/K1x) * (K2x*ux[i]-e*dt/me * Ex[i])
        tEx[i] = Ex[i] + e*NE[i]*dt* (tux[i] + ux[i])/2
        
    ux = copy.deepcopy(tux)
    Ex = copy.deepcopy(tEx)
    scipy.io.savemat('Ex.mat', {'Ex':Ex});
    scipy.io.savemat('ux.mat', {'ux': ux});
    Exex   = map(lambda x: np.exp(complex(0,omega*time))*np.exp(-gamma*x),X)
    Uxex   = map(lambda x: (-e/(me*complex(0,omega)))*x,Exex)
    scipy.io.savemat('Exex.mat', {'Exex':Exex});
    scipy.io.savemat('uxex.mat', {'Uxex':Uxex});

    
        
    temp = 0
    for i in range(N) : 
        temp = temp + np.abs((ux[i]-Uxex[i]))**2
    print 'norm err ux',np.sqrt(dx*temp)

    temp = 0
    for i in range(N) : 
        temp = temp + np.abs((Ex[i]-Exex[i]))**2
    print 'norm err Ex', np.sqrt(dx*temp)    
    
    t = t+1
