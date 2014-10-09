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
N      = 400
# domain ]-L ; L[, (mL=-L)
mL     = -20
pL      =  20
# space step 
dx     = (pL-mL)/N
print "dx", dx
# time step
dt     = 0.5*dx
print "dt",dt
# Number of time steps
Ntime = 100
# constants 
me     = 1#9.11e-31               #electron mass (kg)
e      = -1#1.6e-19              #electron charge (coulombs)
eps0   = 1#8.85e-12               #permeability  (F/m)
nu     = 0                    #friction
B0     = 1#0.95                        #given
c      = 1


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
#wc    = abs(e)*B0/me
#print "omega_c" ,wc
def Ne(x) : 
    if (x < -10) :
        return 0
    elif (x >= -10) and (x <= 10 ):
        return 36*(x+10)/20
    else :
        return 36
    

omega = 3


gamma = 2
#---------------------------------------------------------#

#--------------G for left BC -------------------#
alpha = 1

def g(t) : 
    #return np.exp(complex(0,omega*t))*np.exp(-gamma*mL)*(-gamma+complex(0,alphamL))
    #return 1/2*(1+tanh((t-2)/0.1))*np.exp(-omega*t*complex(0,1))
    return cos(omega*t)
    
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
WP   = map(lambda x   : sqrt(e*e*x/(me*eps0)), NE)
wc   = abs(e)*B0/me
WPy  = map(lambda x   : sqrt(e*e*x/(me*eps0)), NEy)
Wy   = map(lambda x   : x*x + wc*wc, WPy)
print 'wp = ', WP[1]
print 'wc = ', wc
#NE  = map(lambda x : 3.13895515115395 + 3.66439242009228*complex(0,1) ,  X)
#NEy = map(lambda x : 3.13895515115395 + 3.66439242009228*complex(0,1),  X12)
#Ey  = map(lambda x : np.exp(-gamma*x),X12)
#H   = map(lambda x : -(gamma/complex(0,omega))*np.exp(-gamma*x)*np.exp(complex(0,-omega*dt/2)), X)
#uy  = map(lambda x,n : -(1/(e*n))*(eps0*complex(0,omega)-complex(0,gamma*gamma/omega))\
#               *np.exp(-gamma*x), X12,NEy)
#ux   = map(lambda x,n: (1/(B0*e)*(-me/(e*n)*(eps0*omega*omega+gamma*gamma)-e+nu*me/(e*n)*(complex(0,-omega-gamma*gamma/omega))))*np.exp(-gamma*x), X,NE)
#Ex   = map(lambda x,n: ((e*n)/(eps0*complex(0,omega)))*x, ux,NE)

#alpha(-L)
alphamL = 1#(omega*omega)/(c*c)*(1-WP[1]*WP[1]/(omega*omega-wc*wc))


tux = copy.deepcopy(ux)
tuy = copy.deepcopy(uy)
tEx = copy.deepcopy(Ex)
tEy = copy.deepcopy(Ey)
H12 = copy.deepcopy(H)
scipy.io.savemat('H120.mat', {'H120':H12})
T = []
nEx = []
f, axarr = plt.subplots(4,2)
f.subplots_adjust(hspace=.5)
axarr[3,1].axis('off')
plt.tight_layout()
plt.ion()
plt.show()
def plot2(X,T,X2,T2,X3,T3,X4,T4,X5,T5,X6,T6,X7,T7,s,s2,s3):
    axarr[0,0].hold(False)
    axarr[1,0].hold(False)
    axarr[2,0].hold(False)
    axarr[1,1].hold(False)
    axarr[0,1].hold(False)
    axarr[2,1].hold(False)
    axarr[3,0].hold(False)
    
#ylim(-20,20)
    
    axarr[1,0].plot(X, T)
    axarr[1,0].set_title('Ey ')
    axarr[1,1].plot(X2, T2)
    axarr[1,1].set_title('Ex ')
    
    axarr[0,1].plot(X4, T4)
    axarr[0,1].set_title('ux ')
    axarr[0,0].plot(X5, T5)
    axarr[0,0].set_title('uy ')
    axarr[2,1].plot(X6, T6)
    axarr[2,1].set_title('||Ex||')
    axarr[2,0].plot(X3, T3)
    axarr[2,0].set_title('H')
    axarr[3,0].plot(X7, T7)
    axarr[3,0].set_title('Energy')
#axarr.plot(X,np.imag(T),s,label=s2)
#xlabel(r'$x$',fontsize =16)
#ylabel('im',fontsize =16)
#leg = axarr.legend(shadow = True, loc = 3)
    #plt.draw()


#---------------------------------------------------------#
#beta = B0 *dt*e/me
#delta  = e*dt/(me)
#--------------------K1&K2 s def function-------------------#
def Kcoeff(x): 
    #K1     =  1 + nu * dt/2 + dt*dt*x*x/(4)
    #K2     =  1 - nu * dt/2 - dt*dt*x*x/(4)
    #K1x    =  K1 + wc*dt*wc*dt/(4*K1)
    #K2x    =  K2 - wc*dt*wc*dt/(4*K1)
    K1     =  1 + nu * dt/2 + dt*dt*x*e*e/(4*me*eps0)
    K2     =  1 - nu * dt/2 - dt*dt*x*e*e/(4*me*eps0)
    K1x    =  K1 + dt*dt*e*e*B0*B0/(4*K1*me*me)
    K2x    =  K2 -  dt*dt*e*e*B0*B0/(4*K1*me*me)
    return [K1,K2,K1x,K2x]



#-----------------------------------------------------------#
#------------------------time loop--------------------------#
#-----------------------------------------------------------#


t = 1
#coeff alpha (1/alpha_Ey)Ey in left BC
#alphaey = complex(0,alphamL)/2 - 1/dx
#print 'alphaey', alphaey

Et = []
while (t<=Ntime):
    #print "iter : ",t,"time : ", t*dt

    #---------------- H -> H12
    #left BC : Robin BC
    T.append(t*dt)
    H12[0] =  g(t*dt)

    i = 1
    while (i<N):
        H12[i] = H[i] - (dt/dx)*(Ey[i] - Ey[i-1])
        i      = i + 1
    #right BC : Ey[N] = 0
    H12[N] = H[N] +(dt/dx)*(Ey[N-1])

    #------------------- ux -> tux
    for i in range(N):
        [K1,K2,K1x,K2x] = Kcoeff(NEy[i])
        #[K1,K2,K1x,K2x] = Kcoeff(WP[i])

        #tux[i] =  (1/K1x) * (K2x*ux[i] -(wc/B0)*dt*Ex[i] + (((wc*wc/B0)*dt*dt)/(2*K1))*Ey[i]\
        #                         - (wc*wc/B0)*dt*dt*dt/(4*K1)*((H12[i+1] - H12[i])/dx)   \

        #                         - (wc*dt/2)*(K2/K1 + 1)*uy[i])
        tux[i] =  (1/K1x) * (K2x*ux[i] +(e/me)*dt*Ex[i] + (((e*e*B0/(me*me))*dt*dt)/(2*K1))*Ey[i]\
                                 - (e*e*B0/(me*me*eps0))*dt*dt*dt/(4*K1)*((H12[i+1] - H12[i])/dx)   \
                                 + (e*B0*dt/(2*me))*(K2/K1 + 1)*uy[i])
#H[i+1]-H[i] ou H[i]-H[i-1]??
        
    #right BC (Ey(N) = 0, uy(N) = 0)
    #[K1,K2,K1x,K2x] = Kcoeff(WP[N])
    [K1,K2,K1x,K2x] = Kcoeff(NE[N])
    #tux[N] =  (1/K1x)*(K2x*ux[N] -(wc/B0)*dt*Ex[N]\
                           #-(wc*wc/B0)*dt*dt*dt/(4*K1)*(- H12[N])/dx )
                           #+ (((wc*wc/B0)*dt*dt)/(2*K1))*Ey[N-1]\
                           #- (wc*dt/2)*(K2/K1 + 1)*uy[N-1])
    tux[N] =  (1/K1x) * (K2x*ux[N] +(e/me)*dt*Ex[N] \
                             - (e*e*B0/(me*me*eps0))*dt*dt*dt/(4*K1)*(( - H12[N])/dx)  )
    
    #scipy.io.savemat('uy.mat',  {'uy': uy}); 

    #-------------------- u y->tuy

    for i in range(N):
        [K1,K2,K1x,K2x] = Kcoeff(NEy[i])
        #[K1,K2,K1x,K2x] = Kcoeff(WPy[i]);
        #tuy[i] = (1/K1) * (K2 * uy[i] -dt * (wc/B0)* Ey[i] +dt*(dt/2)*(wc/(B0*eps0)) * (H12[i+1] - H12[i])/dx) +wc*dt*(tux[i]+ux[i])/2
        tuy[i] = (1/K1) * (K2 * uy[i] +dt * (e/me)* Ey[i] -dt*(dt/2)*(e/(eps0*me)) * (H12[i+1] - H12[i])/dx) -(e*B0/me)*dt*(tux[i]+ux[i])/2

    

    #------------------- E -> tE
        

    for i in range(N):
        #tEx[i] = Ex[i] + (dt*(WP[i]*WP[i]/wc)*B0)* (tux[i] + ux[i])/(2)
        tEx[i] = Ex[i] - dt*(e*NEy[i]/eps0)* (tux[i] + ux[i])/(2)
        #tEy[i] = Ey[i] - (dt/eps0) * (H12[i+1] - H12[i])/dx  +(dt/2)*(WPy[i]*WPy[i]/wc)*B0*(tuy[i] + uy[i])
        tEy[i] = Ey[i] - (dt/eps0) * (H12[i+1] - H12[i])/dx  -(dt/2)*(e*NEy[i]/eps0)*(tuy[i] + uy[i])
        
    #right BC (nothing to do)
  
    #tEx[N]  = Ex[N] + (dt*(WP[N]*WP[N]/wc)*B0)*(tux[N] + ux[N])/(2) 
    tEx[N] = Ex[N] - dt*(e*NEy[N]/eps0)* (tux[N] + ux[N])/(2)
    
    Ec = 0
    for i in range(N):
        Ec = Ec+NEy[i]*(ux[i]**2+uy[i]**2)

    corr = 0
    
    for i in range(N-2) : 
        corr = corr + (Ey[i+1]-Ey[i])*H12[i+1]
        
    E = (dt/dx)*corr + Ec
    for i in range(N) :
        E = E+ Ey[i]**2+Ex[i]**2 + H12[i]**2 
    


    Et.append(E)
    #print "E",Et

    #----------f^n <- f^n+1
    ux = copy.deepcopy(tux)
    Ex = copy.deepcopy(tEx)
    uy = copy.deepcopy(tuy)
    Ey = copy.deepcopy(tEy)
    H  = copy.deepcopy(H12)
    nEx.append(max(np.absolute(Ex)))
  
    plot2(X12,Ey,X,Ex,X,H,X,ux,T,Et,T,nEx,'b','$E$','E')
    #time = t*dt

    scipy.io.savemat('H.mat',  {'H': H})
    scipy.io.savemat('Ey.mat',  {'Ey': Ey})
    scipy.io.savemat('uy.mat',  {'uy': uy})
    scipy.io.savemat('Ex.mat',  {'Ex': Ex})
    scipy.io.savemat('ux.mat',  {'ux': ux})
    scipy.io.savemat('X.mat', {'X': X})
    scipy.io.savemat('X2.mat', {'X2': X12})
    scipy.io.savemat('T',{'T':T})
    scipy.io.savemat('nEx',{'nEx':nEx})
    t = t+1           



