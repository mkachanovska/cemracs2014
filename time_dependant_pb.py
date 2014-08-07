#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import os
from time import time
import numpy as np
from pylab import *
from scipy.special import airy


#---------------CONSTANTS--------------------#
# mesh size
N      = 124
# domain ]-L ; H[, (mL=-L)
mL     = -20
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
nu     = 1                     #friction
B0     = 1                      #given
c      = 1#3e8 


#gamma/beta/nudelta
gamma  = -e*dt / (me*c*c*eps0)
beta   = -abs(e)*B0*dt/(c*me)
nud    = nu*dt/c
print 'gamma',gamma
print 'beta',beta
print 'nudelta', nud
#alpha(-L)
alphamL = 2

#coeff alpha (1/alpha_Ey)Ey in left BC
alphaey = complex(0,alphamL)/2 - 1/dx
print 'alphaey', alphaey
#-------------airy functions------------------#
#return A
def airy_one(x):
    [A,B,C,D] = airy(x)
    return A
#return A'
def airy_two(x):
    [A,B,C,D] = airy(x)
    return B
#-----------------------------------------------#

#--------------G for left BC -------------------#
def g(t) : 
    return (airy_two(mL)+complex(0,alphamL)*airy_one(mL))*np.exp(complex(0,t))*eps0*c

#-----------------------------------------------#
#------table and functions initialisation-------#
#-----------------------------------------------#
X    = map(lambda i: mL + i*dx, range(N+2))
X12  = map(lambda i: mL+0.5*dx + i*dx, range(N+1))
#Ex   = map(lambda x: np.exp(complex(0,k*x)),X)
Ex   = map(lambda x: 0, X)
#ux   = map(lambda x: np.exp(complex(0,k*x)),X)
ux   = map(lambda x: 0, X)
#H    = map(lambda x: np.exp(complex(0,k*x)),X)
H    = map(lambda x: 0, X)
#Ey   = map(lambda x: np.exp(complex(0,k*x)),X12)
Ey   = map(lambda x: airy_one(x),X12)
AiryVec = map(lambda x: airy_one(x),X12)
#uy   = map(lambda x: np.exp(complex(0,k*x)),X12)
uy    = map(lambda x: 0, X)

H12 = H
tux = ux
tuy = uy
tEx = Ex
tEy = Ey
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

NE    = map(lambda x: Ne(x), X)

#---------------------------------------------------------#

#----------------------- omega_p -------------------------#
def wp(x) : 
    return np.sqrt(e*e*Ne(x)/(me*eps0))


#---------------------------------------------------------#

#---------plot functions definition (real, imag)----------#
def plot_real(X,T,s,s2,s3):
    
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    #ax.set_xlim(-np.pi-0.01,np.pi)
    plot(X,real(T),s,label=s2)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    xlim([-20,10])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()
def plot_imag(X,T,s,s2,s3):
    fig = plt.figure() 
    ax = fig.add_subplot(111)   
    plot(X,np.imag(T),s,label=s2)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    xlim([-20,10])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()

    

#--------------------K1&K2 s def function-------------------#
def Kcoeff(x): 
    K1     =  1 + nud/2. - (gamma*dt)*e*x/4
    K2     =  1 - nud/2. + (gamma*dt)*e*x/4
    K1x    =  K1 + beta*beta/(4*K1)
    K2x    =  K2 - beta*beta/(4*K1)
    return [K1,K2,K1x,K2x]
#-----------------------------------------------------------#
#------------------------time loop--------------------------#
#-----------------------------------------------------------#

for t in range(Ntime):
    

    #------------------- ux -> tux
    for i in range(N+1):
        
        [K1,K2,K1x,K2x] = Kcoeff(NE[i])
        tux[i] =  (1/K1x) * (K2x*ux[i] + gamma*Ex[i] + ((beta*gamma)/(2*K1))*Ey[i]\
                               - beta*gamma*dt/(4*K1)*((H12[i] - H12[i-1])/dx)   \
                               + (beta/2)*(K2/K1 + 1)*uy[i])
        
    #left BC (if i = N) (H(i+1) = 0) ????
    [K1,K2,K1x,K2x] = Kcoeff(NE[N+1])
    tux[N+1] =  (1/K1x)*(K2x*ux[N+1] + gamma*Ex[N+1] + ((beta*gamma)/(2*K1))*Ey[N]\
                              - beta*gamma*dt/(4*K1)*(H12[N+1] - H12[N])/dx   \
                              + (beta/2)*(K2/K1 + 1)*uy[N])

    #-------------------- uy->tuy
    for i in range(N):
        [K1,K2,K1x,K2x] = Kcoeff(NE[i])
        tuy[i] =  (1/K1)*(K2*uy[i] + gamma*Ey[i] - gamma*dt/2 *(H12[i+1] - H12[i])/dx\
                              - beta*(tux[i] + ux[i])/2)
    #left BC (if i = N)(H(i+1) = 0) ????
    [K1,K2,K1x,K2x] = Kcoeff(NE[N+1])
    tuy[N] =  (1/K1)*(K2*uy[N] + gamma*Ey[N] - gamma*dt/2 *(H12[N+1] - H12[N])/dx\
                          - beta*(tux[N] + ux[N])/2)
        
    #------------------- E -> tE
    for i in range(N):
        tEx[i] = Ex[i] + e*NE[i]*dt* (tux[i] + ux[i])/2
        tEy[i] = Ey[i] - (dt) * (H12[i+1] - H12[i])/dx \
            + (dt*e*NE[i]/2)*(tuy[i] + uy[i])
        
   #left BC (if i = N)(H(i+1) = 0) ????
    tEx[N]    = Ex[N] - (dt*e*NE[N])*(tux[N] + ux[N])/2
    tEx[N+1]  = Ex[N+1] - (dt*e*NE[N+1])*(tux[N+1] + ux[N+1])/2
    tEy[N]    = Ey[N] - (dt) * (H12[N+1]- H12[N])/dx \
        - (dt*e*NE[N]/2)*(tuy[N] + uy[N])
    #---------------- H -> H12
    H12[0] = H[0] -  dt/dx * (Ey[0] \
                                  - (1/alphaey) * (g(t*dt+dt/2) -  complex(0,alphamL)/2 * Ey[0] \
                                                       - (1/dx)* Ey[0]))
    i = 1
    while (i<N+1):
        im12   = i-1
        i12    = i+1
        H12[i] = H[i] - dt/dx*(Ey[i-1] - Ey[i])
        i      = i +1
    H12[N+1] = H[N+1] 


    #----------f^n <- f^n+1
    ux = tux
    Ex = tEx
    uy = tuy
    Ey = tEy
    H  = H12

    print "max(ux)" (max(ux))
    print "max(Ex)" (max(Ex))
    print "max(uy)" (max(uy))
    print "max(Ey)" (max(Ey))
    print "max(Hz)" (max(H))


#-------plots------#
plot_real(X,ux,'b-','ux','ux')
plot_imag(X,ux,'b-','ux','ux')
