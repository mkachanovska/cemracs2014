#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import os
import time
import numpy as np
from pylab import *
import copy
import sub
#---------------CONSTANTS--------------------#
# mesh size
N      = 2048
# domain ]-L ; L[, (mL=-L)
mL     = -20
pL      = 10
# space step 
dx     = (pL-mL)/(N-1)
print "dx", dx
# time step
dt     = 0.5*dx
print "dt",dt
# Number of time steps
Ntime = 2e6
# constants 
me     = 1#9.11e-31               #electron mass (kg)
e      = -1#1.6e-19              #electron charge (coulombs)
eps0   = 1#8.85e-12               #permeability  (F/m)
nu     = 1e4                  #friction
B0     = 0#0.95                        #given
c      = 1


#-------------------------- NE --------------------------#

def Ne(x) : 
	return 0.2;
#    if (x < -10) :
#        return 0
#    elif (x >= -10) and (x <= 10 ):
#        return 36*(x+10)/20
#    else :
#        return 36
    

omega = 1


gamma = 2
#---------------------------------------------------------#

#--------------G for left BC -------------------#
alpha = 1

def g(t) : 
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


alphamL = 1

#f, axarr = plt.subplots(4,2)
#f.subplots_adjust(hspace=.5)
#axarr[3,1].axis('off')
#plt.tight_layout()
#plt.ion()
#plt.show()
#def plot2(X,T,X2,T2,X3,T3,X4,T4,X5,T5,X6,T6,X7,T7,s,s2,s3):
#    axarr[0,0].hold(False)
#    axarr[1,0].hold(False)
#    axarr[2,0].hold(False)
#    axarr[1,1].hold(False)
#    axarr[0,1].hold(False)
#    axarr[2,1].hold(False)
#    axarr[3,0].hold(False)
#    
#ylim(-20,20)
#    
#    axarr[1,0].plot(X, T)
#    axarr[1,0].set_title('Ey ')
#    axarr[1,1].plot(X2, T2)
#    axarr[1,1].set_title('Ex ')
#    
#    axarr[0,1].plot(X4, T4)
#    axarr[0,1].set_title('ux ')
#    axarr[0,0].plot(X5, T5)
#    axarr[0,0].set_title('uy ')
#    axarr[2,1].plot(X6, T6)
#    axarr[2,1].set_title('||Ex||')
#    axarr[2,0].plot(X3, T3)
#    axarr[2,0].set_title('H')
#    axarr[3,0].plot(X7, T7)
#    axarr[3,0].set_title('Energy')




#-----------------------------------------------------------#
#------------------------time loop--------------------------#
#-----------------------------------------------------------#




sub.tdp.tdp_sub(dx,dt,Ntime,me,e,eps0,nu,B0,omega,NE,NEy,ux,H,Ex,Ey,uy,X,X12)


#plot2(X12,Ey,X,Ex,X,H,X,ux,T,Et,T,nEx,'b','$E$','E')




