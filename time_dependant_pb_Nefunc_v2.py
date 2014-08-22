#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import os
import time
import numpy as np
from pylab import *
import scipy.io
from scipy.special import airy
import copy

#---------------CONSTANTS--------------------#
# mesh size
N      = 500
# domain ]-L ; H[, (mL=-L)
mL     = -1
H      =  200
# space step 
dx     = (H-mL)/N
print "dx", dx
# time step
dt     = 0.5*dx
print "dt",dt
# Number of time steps
Ntime = 5000
# constants 
me     = 1#9.11e-31               #electron mass (kg)
e      = -1#1.6e-19              #electron charge (coulombs)
eps0   = 1#8.85e-12               #permeability  (F/m)
nu     = 0                       #friction
B0     = 1 #0.95                        #given
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
wc    = abs(e)*B0/me
print "omega_c" ,wc
def Ne(x) : 
    return 1#(2*B0)/(1+exp(-x)) ;
    

omega = 0
wp = 1

gamma = 2
#---------------------------------------------------------#

#--------------G for left BC -------------------#
alpha = 1

def g(t) : 
    #return np.exp(complex(0,omega*t))*np.exp(-gamma*mL)*(-gamma+complex(0,alphamL))
    return -1/2*(1+tanh((t-2)/0.1))*np.exp(omega*t*complex(0,1)) 

print g(0)
#-----------------------------------------------#
#------table and functions initialisation-------#
#-----------------------------------------------#
X    = map(lambda i : mL + i*dx, range(N+1))
NE   = map(lambda x : Ne(x) ,  X)
X12  = map(lambda i : mL+0.5*dx + i*dx, range(N))
NEy  = map(lambda x : Ne(x),X12)
Ey   = map(lambda x : 0, X12)
H    = map(lambda x : 0, X)
uy   = map(lambda x : 0, X12)
ux   = map(lambda x : 0, X) 
Ex   = map(lambda x : 0, X)

#NE  = map(lambda x : 3.13895515115395 + 3.66439242009228*complex(0,1) ,  X)
#NEy = map(lambda x : 3.13895515115395 + 3.66439242009228*complex(0,1),  X12)
#Ey  = map(lambda x : np.exp(-gamma*x),X12)
#H   = map(lambda x : -(gamma/complex(0,omega))*np.exp(-gamma*x)*np.exp(complex(0,-omega*dt/2)), X)
#uy  = map(lambda x,n : -(1/(e*n))*(eps0*complex(0,omega)-complex(0,gamma*gamma/omega))\
#               *np.exp(-gamma*x), X12,NEy)
#ux   = map(lambda x,n: (1/(B0*e)*(-me/(e*n)*(eps0*omega*omega+gamma*gamma)-e+nu*me/(e*n)*(complex(0,-omega-gamma*gamma/omega))))*np.exp(-gamma*x), X,NE)
#Ex   = map(lambda x,n: ((e*n)/(eps0*complex(0,omega)))*x, ux,NE)


tux = copy.deepcopy(ux)
tuy = copy.deepcopy(uy)
tEx = copy.deepcopy(Ex)
tEy = copy.deepcopy(Ey)
H12 = copy.deepcopy(H)
scipy.io.savemat('H120.mat', {'H120':H12})



#---------------------------------------------------------#
beta = B0 *dt/me
delta  = e*dt/(me)
#--------------------K1&K2 s def function-------------------#
def Kcoeff(x): 
    K1     =  1 + nu * dt/2 + dt*dt*x*e*e/(4*me*eps0)
    K2     =  1 - nu * dt/2 - dt*dt*x*e*e/(4*me*eps0)
    K1x    =  K1 + beta*beta/(4*K1)
    K2x    =  K2 - beta*beta/(4*K1)
    
    return [K1,K2,K1x,K2x]



#-----------------------------------------------------------#
#------------------------time loop--------------------------#
#-----------------------------------------------------------#
fig = plt.figure() 
ax = fig.add_subplot(111)

plt.ion()
plt.show()
def plot_real(X,T,s,s2,s3):
    #ax.set_xlim(mL,H)
    clf()
    plot(X,real(T),s,label=s2)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    leg = ax.legend(shadow = True, loc = 3)
    plt.draw()
    time.sleep(0.05)
def plot_imag(X,T,s,s2,s3):
    #ax.set_xlim(mL,H)

    clf()
    plot(X,np.imag(T),s,label=s2)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    leg = ax.legend(shadow = True, loc = 3)
    plt.draw()
    time.sleep(0.05)



t = 1
#coeff alpha (1/alpha_Ey)Ey in left BC
alphaey = complex(0,alphamL)/2 - 1/dx
print 'alphaey', alphaey


while (t<=Ntime):
    print "iter : ",t,"time : ", t*dt

    #---------------- H -> H12
    H12[0] = H[0] -  dt/dx* (Ey[0] - 1/alphaey * (Ey[0] *(- 1/dx - complex(0,alphamL)/2)+g((t-1)*dt)))

    i = 1
    while (i<N):
        H12[i] = H[i] - dt/dx*(Ey[i] - Ey[i-1])
        i      = i + 1
    H12[N] = H[N]

    #------------------- ux -> tux
    for i in range(N):
        
        [K1,K2,K1x,K2x] = Kcoeff(NE[i])
        tux[i] =  (1/K1x) * (K2x*ux[i] + delta*Ex[i] + ((beta*delta)/(2*K1))*Ey[i]\
                                 - beta*delta*dt/(4*K1)*((H12[i] - H12[i-1])/dx)   \
                                 +(beta/2)*(K2/K1 + 1)*uy[i])
        
    #left BC (if i = N) (H(i+1) = 0) ????
    [K1,K2,K1x,K2x] = Kcoeff(NE[N])
    tux[N] =  (1/K1x)*(K2x*ux[N] + delta*Ex[N]\
                            -beta*delta*dt/(4*K1)*(H12[N] - H12[N-1])/dx\
                            +((beta*delta)/(2*K1))*Ey[N-1]\
                            +(beta/2)*(K2/K1 + 1)*uy[N-1])

    scipy.io.savemat('uy.mat',  {'uy': uy});
 
    #plot_real(X,tux,'b','$u_x$','u_x')
    #-------------------- u y->tuy

    for i in range(N):
        [K1,K2,K1x,K2x] = Kcoeff(NEy[i]);
        tuy[i] = (1/K1) * (K2 * uy[i] + delta * Ey[i] -dt*dt*e/(2*me*eps0) * (H12[i+1] - H12[i])/dx)- beta*(tux[i]+ux[i])/2
    #plot_real(X12,tuy,'b','$u_y$','u_y')
    
    #------------------- E -> tE
                        

    for i in range(N):
        tEx[i] = Ex[i] - e*NE[i]*dt* (tux[i] + ux[i])/(2*eps0)
        tEy[i] = Ey[i] - (dt/eps0) * (H12[i+1] - H12[i])/dx  -(dt*e*NEy[i]/(2*eps0))*(tuy[i] + uy[i])
    
    #left BC (if i = N)(H(i+1) = 0) ????
  
    tEx[N]  = Ex[N] - (dt*e*NE[N])*(tux[N] + ux[N])/(2*eps0)

    #----------f^n <- f^n+1
    ux = copy.deepcopy(tux)
    Ex = copy.deepcopy(tEx)
    uy = copy.deepcopy(tuy)
    Ey = copy.deepcopy(tEy)
    H  = copy.deepcopy(H12)

    

    scipy.io.savemat('H.mat',  {'H': H})
    scipy.io.savemat('results/Ey'+str(t)+'.mat',  {'Ey': Ey})
    scipy.io.savemat('uy.mat',  {'uy': uy})
    scipy.io.savemat('Ex.mat',  {'Ex': Ex})
    scipy.io.savemat('ux.mat',  {'ux': ux})

    t = t+1           



