#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import os
from time import time
import numpy as np
from pylab import *
from scipy.special import airy
#mesh size
N = 256
#domain ]-L ; H[, (mL=-L)
mL = -20
H =  10
#space step 
dx = (H-mL)/N
#time step
dt = 0.5*dx
#final time 
tfinal = 10
#constants 
me = 9.11e-31#electron mass (kg)
Ne = 1e19 #density of e- f(x)
e =  -1.6e-19 #electron charge (coulombs)
eps0 = 8.85e-12 #permeability  (F/m)
nu = 0#friction
B0 = 1 #given
c = 3e8 

gamma = e*dt / (me*c*c*eps0)
beta = abs(e)*B0*dt/(c*me)
nud = nu*dt/c

alphamL = 2
alphaey = complex(0,alphamL)/2 - 1/dx
#airy function return
def airy_one(x):
    [A,B,C,D] = airy(x)
    return A
def airy_two(x):
    [A,B,C,D] = airy(x)
    return B

#initialisation
def g(t) : 
    return (airy_two(mL)+complex(0,alphamL)*airy_one(mL))*np.exp(complex(0,t))*eps0*c
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
wc= abs(e*B0)/me
kappa = wc/(c*(c*c-wc*wc))
psi = 1/(c*c-wc*wc)
def Ne(x) : 
    a = (kappa*kappa- psi*psi)*e*e*e*e/(me*me*eps0*eps0)
    b =  (x*psi+2*psi)*(e/(me*eps0))
    c = x-1
    delta = (x*psi+2*psi-psi*psi)**2 - 4*(x-1)*(kappa*kappa*e*e*e*e/(me*me*eps0*eps0))
    print delta 
    if (delta == 0):
        return -b/(2*a)
    elif (delta >0):
        return (-b + np.sqrt(delta))/(2*a)
    else : 
        return (-b +complex(0,1)*np.sqrt(-delta))/(2*a)

def wp(x) : 
    return np.sqrt(e*e*Ne(x)/(me*eps0))
    

def plot_real():
    
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    #ax.set_xlim(-np.pi-0.01,np.pi)
    plot(X12,real(Ey),'b-',label=r"Re(E_y)")
    plot(X12,real(AiryVec),"r--", label=r"Re(airy)")
    xlabel(r'$x$',fontsize =16)
    ylabel(r'$E$',fontsize =16)
    xlim([-20,10])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()
def plot_imag():
    fig = plt.figure() 
    ax = fig.add_subplot(111)   
    plot(X12,np.imag(Ey),'c-',label=r"Im(E_y)")
    plot(X12,np.imag(AiryVec),"m--", label=r"Im(airy)")
    xlabel(r'$x$',fontsize =16)
    ylabel(r'$E$',fontsize =16)
    xlim([-20,10])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()

#time loop
H12 = H
tux = ux
tuy = uy
tEx = Ex
tEy = Ey
for t in range(tfinal):
    

    # ux -> tux
    for i in range(N+1):
        
        
        K1 =  1 + nud/2. + (gamma*dt)*e*NE[i]/4
        K2 =  1 - nud/2. - (gamma*dt)*e*NE[i]/4
        K1x = K1 + beta*beta/(4*K1)
        K2x = K2-beta*beta/(4*K1)
        tux[i] = (1/K1x)*(K2x*ux[i]+gamma*Ex[i]+((beta*gamma)/(2*K1))*Ey[i]\
                              - beta*gamma*dt/(4*K1)*((H12[i]-H12[i-1])/dx)   \
                              + (beta/2)*(K2/K1+1)*uy[i])
        
    #left BC (if i = N) (H(i+1) = 0) ????
    K1 =  1 + nud/2. + (gamma*dt)*e*NE[N+1]/4
    K2 =  1 - nud/2. - (gamma*dt)*e*NE[N+1]/4
    K1x = K1 + beta*beta/(4*K1)
    K2x = K2-beta*beta/(4*K1)
    tux[N+1] =  (1/K1x)*(K2x*ux[N+1]+gamma*Ex[N+1]+((beta*gamma)/(2*K1))*Ey[N]\
                              - beta*gamma*dt/(4*K1)*(H12[N+1]-H12[N])/dx   \
                              + (beta/2)*(K2/K1+1)*uy[N])

    #uy->tuy
    for i in range(N):
        K1 =  1 + nud/2. + (gamma*dt)*e*NE[i]/4
        K2 =  1 - nud/2. - (gamma*dt)*e*NE[i]/4
        K1x = K1 + beta*beta/(4*K1)
        K2x = K2-beta*beta/(4*K1)
        tuy[i] = (1/K1)*(K2*uy[i]+ gamma*Ey[i] - gamma*dt/2 *(H12[i+1]-H12[i])/dx\
                             -beta*(tux[i]+ux[i])/2)
    #left BC (if i = N)(H(i+1) = 0) ????
    K1 =  1 + nud/2. + (gamma*dt)*e*NE[N+1]/4
    K2 =  1 - nud/2. - (gamma*dt)*e*NE[N+1]/4
    K1x = K1 + beta*beta/(4*K1)
    K2x = K2-beta*beta/(4*K1)
    tuy[N] = (1/K1)*(K2*uy[N]+ gamma*Ey[N] - gamma*dt/2 *(H12[N+1]-H12[N])/dx\
                             -beta*(tux[N]+ux[N])/2)
        
    #E -> tE
    for i in range(N):
        tEx[i] = Ex[i] -e*NE[i]*dt* (tux[i]+ux[i])/2
        tEy[i] = Ey[i] - (dt) * (H12[i+1] - H12[i])/dx \
            - (dt*e*NE[i]/2)*(tuy[i]+uy[i])
        
   #left BC (if i = N)(H(i+1) = 0) ????
    tEx[N] = Ex[N] - (dt*e*NE[N])*(tux[N]+ux[N])/2
    tEx[N+1] = Ex[N+1] - (dt*e*NE[N+1])*(tux[N+1]+ux[N+1])/2
    tEy[N] = Ey[N] - (dt) * (H12[N+1]- H12[N])/dx \
        - (dt*e*NE[N]/2)*(tuy[N]+uy[N])
    # H advance in dt/2
    H12[0] = H[0] -  dt/dx * (Ey[0] \
                                  - (1/alphaey) * (g(t*dt+dt/2) -  complex(0,alphamL)/2 * Ey[0] \
                                                       - (1/dx)* Ey[0]))
    i = 1
    while (i<N+1):
        im12 = i-1
        i12 = i+1
        H12[i] = H[i] - dt/dx*(Ey[i-1]-Ey[i])
        i = i +1
    H12[N+1] = H[N+1] 
    ux = tux
    Ex = tEx
    uy = tuy
    Ey = tEy
    H  = H12

plot_real()
plot_imag()
