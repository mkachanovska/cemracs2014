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
N      = 500
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


#gamma/beta/nudelta
#gamma  = -e*dt / (me*c*c*eps0)

#beta   = -abs(e)*B0*dt/(c*me)
#nud    = nu*dt/c
#print 'gamma',gamma
#print 'beta',beta
#print 'nudelta', nud
#alpha(-L)
alphamL = 2

#coeff alpha (1/alpha_Ey)Ey in left BC
alphaey = complex(0,alphamL)/2 - 1/dx
print 'alphaey', alphaey
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
    return np.exp(complex(0,omega*t))*np.exp(-gamma*mL)*(-gamma+complex(0,alphamL));
#(1/dt)*np.exp(complex(0,omega*t)-gamma*mL)*(1-np.exp(-gamma*dx))+complex(0,alphamL)*np.exp(-gamma*mL+complex(0,omega*t)) 
#(airy_two(mL)+complex(0,alphamL)*airy_one(mL))*np.exp(complex(0,t))*eps0*c
#-----------------------------------------------#
#------table and functions initialisation-------#
#-----------------------------------------------#
X    = map(lambda i: mL + i*dx, range(N+1))
NE    = map(lambda x: 3*me/(e*e),  X)
X12  = map(lambda i: mL+0.5*dx + i*dx, range(N))
#Ex   = map(lambda x: np.exp(complex(0,k*x)),X)
Ex   = map(lambda x: 0, X)
#ux   = map(lambda x: np.exp(complex(0,k*x)),X)
ux   = map(lambda x: 0, X)
#H    = map(lambda x: np.exp(complex(0,k*x)),X)
H    = map(lambda x: (gamma/complex(0,omega))*np.exp(-gamma*x), X)
Ey   = map(lambda x: np.exp(-gamma*x),X12)
#Ey   = map(lambda x: 0,X12)
#AiryVec = map(lambda x: airy_one(x),X12)
#uy   = map(lambda x: np.exp(complex(0,k*x)),X12)
uy    = map(lambda x: -e/(me*complex(0,omega))*np.exp(-gamma*x), X12)

#H12 = map(lambda x: -np.exp(-alpha*dt/2)*np.exp(-alpha*x),X)

tux = copy.deepcopy(ux)
tuy = copy.deepcopy(uy)
tEx = copy.deepcopy(Ex)
tEy = copy.deepcopy(Ey)
H12 = copy.deepcopy(H)

#----------------------- omega_p -------------------------#
#def wp(x) : 
#    return np.sqrt(e*e*Ne(x)/(me*eps0))


#---------------------------------------------------------#

#---------plot functions definition (real, imag)----------#
def plot_real(X,T,s,s2,s3):
    
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    #ax.set_xlim(-np.pi-0.01,np.pi)
    plot(X,real(T),s,label=s2)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    xlim([-1,50])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()
def plot_imag(X,T,s,s2,s3):
    fig = plt.figure() 
    ax = fig.add_subplot(111)   
    plot(X,np.imag(T),s,label=s2)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    xlim([-1,50])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()

def plot_real2(X,T,U,s,s5,s2,s3,s4):
    
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    #ax.set_xlim(-np.pi-0.01,np.pi)
    plot(X,real(T),s,label=s2)
    plot(X,real(U),s5,label=s4)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    xlim([-1,50])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()

def plot_imag2(X,T,U,s,s5,s2,s3,s4):

    fig = plt.figure() 
    ax = fig.add_subplot(111)   
    plot(X,np.imag(T),s,label=s2)
    plot(X,np.imag(U),s5,label=s4)
    xlabel(r'$x$',fontsize =16)
    ylabel(s3,fontsize =16)
    xlim([-1,50])
    leg = ax.legend(shadow = True, loc = 3)
    plt.show()
    clf()

#--------------------K1&K2 s def function-------------------#
def Kcoeff(x): 
    #K1     =  1 + nud/2. - (gamma*dt)*e*x/4
    K1     =  1 + nu * dt/2 + dt*dt*x*e*e/(4*me)
    #K2     =  1 - nud/2. + (gamma*dt)*e*x/4
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

while (t<=Ntime):
    #---------------- H -> H12
    H12[0] = H12[0] -  dt/dx * (tEy[0] \
                                  - (1/alphaey) * (g(t*dt-dt) -  complex(0,alphamL)/2 * tEy[0] \
                                                       - (1/dx)* tEy[0]))

    #H12[0] = H12[0] - dt/dx *(tEy[0] - np.exp(-alpha*(mL-dx/2))*np.exp(-alpha*t*dt))
    i = 1
    while (i<N):
        H12[i] = H12[i] - dt/dx*(tEy[i] - tEy[i-1])
        i      = i +1
    H12[N] = H12[N]
    scipy.io.savemat('H12.mat', {'H12': H12});

    #------------------- ux -> tux
    for i in range(N-1):
        
       # [K1,K2,K1x,K2x] = Kcoeff(NE[i])
        tux[i] =  0 #(1/K1x) * (K2x*ux[i] + gamma*Ex[i] + ((beta*gamma)/(2*K1))*Ey[i]\
                    #           - beta*gamma*dt/(4*K1)*((H12[i] - H12[i-1])/dx)   \
                    #           + (beta/2)*(K2/K1 + 1)*uy[i])
        
    #left BC (if i = N) (H(i+1) = 0) ????
  #  [K1,K2,K1x,K2x] = Kcoeff(NE[N+1])
    tux[N] =  0#(1/K1x)*(K2x*ux[N+1] + gamma*Ex[N+1] + ((beta*gamma)/(2*K1))*Ey[N]\
                 #             - beta*gamma*dt/(4*K1)*(H12[N+1] - H12[N])/dx   \
                 #             + (beta/2)*(K2/K1 + 1)*uy[N])
    scipy.io.savemat('uy.mat',  {'uy': uy}); 

    #-------------------- u y->tuy
    for i in range(N-2):
        [K1,K2,K1x,K2x] = Kcoeff(NE[i])
        #tuy[i] = 0# (1/K1)*(K2*uy[i] + gamma*Ey[i] - gamma*dt/2 *(H12[i+1] - H12[i])/dx\
                  #            - beta*(tux[i] + ux[i])/2)
        tuy[i] = (1/K1) * (K2 * uy[i] - dt * e / me * Ey[i] + dt*dt*e/(2*me) * (H12[i+1] - H12[i])/dx)
    #left BC (if i = N)(H(i+1) = 0) ????
  #  [K1,K2,K1x,K2x] = Kcoeff(NE[N+1])
    #tuy[N-1] =  0#(1/K1)*(K2*uy[N] + gamma*Ey[N] - gamma*dt/2 *(H12[N+1] - H12[N])/dx\
               #           - beta*(tux[N] + ux[N])/2)
    tuy[N-1] = (1/K1) * (K2 * uy[N-1] - dt * e / me * Ey[N-1] + dt*dt*e/(2*me) * (H12[N] - H12[N-1])/dx)
    #------------------- E -> tE
    scipy.io.savemat('H12.mat',  {'H12': H12});
    scipy.io.savemat('Eyo.mat',  {'Eyo': Ey});
    scipy.io.savemat('tuy.mat',  {'tuy': tuy});                          

    for i in range(N-1):
        tEx[i] = Ex[i] + e*NE[i]*dt* (tux[i] + ux[i])/2

        tEy[i] = Ey[i] - (dt) * (H12[i+1] - H12[i])/dx# \
           # + (dt*e*NE[i]/2)*(tuy[i] + uy[i])

   #left BC (if i = N)(H(i+1) = 0) ????
  
    tEx[N-1]    = Ex[N-1] - (dt*e*NE[N-1])*(tux[N-1] + ux[N-1])/2
    tEx[N]  = Ex[N] - (dt*e*NE[N])*(tux[N] + ux[N])/2
    #tEy[N-1]    = Ey[N-1] - (dt) * (H12[N]- H12[N-1])/dx #\
        #- (dt*e*NE[N]/2)*(tuy[N] + uy[N])
    scipy.io.savemat('Ey.mat',  {'arr': tEy});
     
    

    #----------f^n <- f^n+1
    ux = copy.deepcopy(ux)
    Ex = copy.deepcopy(tEx)
    uy = copy.deepcopy(tuy)
    Ey = copy.deepcopy(tEy)
    H  = copy.deepcopy(H12)
    time = t*dt
    Eyex   = map(lambda x: np.exp(complex(0,omega*time))*np.exp(-gamma*x),X12)
    Hex    = map(lambda x: (gamma/complex(0,omega))*np.exp(-gamma*x)*np.exp(complex(0,omega*(time+dt/2))),X)
    Uyex   = map(lambda x: -(e/me)*np.exp(complex(0,omega*time)-gamma*x),X12)
    temp = 0
    for i in range(N):
        temp = temp + np.abs((H[i]-Hex[i]))
    print time, 'norm err H', np.sqrt(dx*temp)
    temp = 0
    for i in range(N) : 
        temp = temp + np.abs((uy[i]-Uyex[i]))
    print 'norm err uy',np.sqrt(dx*temp)
    
    temp = 0
    for i in range(N) : 
        temp = temp + np.abs((Ey[i]-Eyex[i]))
    print 'norm err Ey', np.sqrt(dx*temp)                         
    #if t==Ntime:
    #    plot_real2(X12,Ey,Eyex,'b-','r-','Ey','re',"Eyex")
    #    plot_imag2(X12,Ey,Eyex,'b-','r-','Ey','im',"Eyex")
    #    plot_real2(X,H ,Hex,'b-','r-','H','re',"Hex")
    #    plot_imag2(X,H ,Hex,'b-','r-','H','im',"Hex")
    t = t+1           
    


#-------plots------#
#plot_real(X,ux,'b-','ux','Re')
#plot_imag(X,ux,'b-','ux','Im')
#plot_real(X12,uy,'b-','uy','Re')
#plot_imag(X12,uy,'b-','uy','Im')
#plot_real(X,Ex,'b-','Ex','Re')
#plot_real(X12,Ey,'b-','Ey','Re')
#plot_imag(X12,Ey,'b-','Ey','Im')
#plot_imag(X,Ex,'b-','Ex','Im')
#plot_real(X,H,'b-','H','Re')

#plot_real2(X12,Ey,Eyex,'b-','r-','Ey','re',"Eyex")
#plot_imag2(X12,Ey,Eyex,'b-','r-','Ey','im',"Eyex")
#plot_real2(X,H ,Hex,'b-','r-','H','re',"Hex")
#plot_imag2(X,H ,Hex,'b-','r-','H','im',"Hex")
