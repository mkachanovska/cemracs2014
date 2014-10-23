#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from pylab import *
import copy
import sub
Niter = 10000

f, axarr = plt.subplots(3,2)
f.subplots_adjust(hspace=.5)
plt.tight_layout()



def plot2(X,ux,Ex,H,X12,uy,Ey,T,Et):
    axarr[0,0].hold(False)
    axarr[1,0].hold(False)
    axarr[2,0].hold(False)
    axarr[1,1].hold(False)
    axarr[0,1].hold(False)
    axarr[2,1].hold(False)
    
    
    ylim(-20,20)
    
    axarr[1,0].plot(X12, Ey)
    axarr[1,0].set_title('Ey ')

    
    axarr[1,1].plot(X, Ex)
    
    axarr[1,1].set_title('Ex ')
    axarr[0,1].plot(X, ux)
    
    axarr[0,1].set_title('ux ')
    axarr[0,0].plot(X12,uy )
    axarr[0,0].set_title('uy ')
    
    axarr[2,0].plot(X, H)
    axarr[2,0].set_title('H')
    
    axarr[2,1].plot(T, Et)
    
    axarr[2,1].set_title('Energy')


iter = 0
while (iter <Niter):
    print 'iter',iter
    mat = np.loadtxt('/UMA/tmp/maryna/password//ux%4.4d.data'%iter)
    X = mat[:,0]
    ux = mat[:,1]
    print 'X', np.isnan(np.sum(X))
    print 'ux', np.isnan(np.sum(ux))
    
    mat = np.loadtxt('/UMA/tmp/maryna/password//Ex%4.4d.data'%iter)
    Ex = mat[:,1]
    print 'Ex', np.isnan(np.sum(Ex))
    
    mat = np.loadtxt('/UMA/tmp/maryna/password//H%4.4d.data'%iter)
    H = mat[:,1]
    print 'H', np.isnan(np.sum(H))

    
    mat = np.loadtxt('/UMA/tmp/maryna/password//uy%4.4d.data'%iter)
    X12 = mat[:,0]
    uy = mat[:,1]
    print 'X12', np.isnan(np.sum(X12))
    print 'uy', np.isnan(np.sum(uy))
    mat = np.loadtxt('/UMA/tmp/maryna/password//Ey%4.4d.data'%iter)
    Ey = mat[:,1]
    print 'Ey',np.isnan(np.sum(Ey))

    mat = np.loadtxt('/UMA/tmp/maryna/password//ET.data')
    T = mat[:,0]
    Et = mat[:,1]
    print "T", np.isnan(np.sum(T))
    print 'Et', np.isnan(np.sum(Ey))
    plot(X,ux)
    plot2(X,ux,Ex,H,X12,uy,Ey,T,Et)
    savefig('im/res%4.4d.png'%iter)
   #quit()
    iter = iter+1

#os.system('avconv -r 10 -i im/res%4d.png test.mp4')
