#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
'''
This is a demo script to see what changing the
parameters of the DTD function-- the unimodal
model that I currently use
L. Strolger
'''


if __name__=='__main__':

    xx = arange(0,10,0.1)

    ax1 = subplot(131)
    ax2 = subplot(132)
    ax3 = subplot(133)

    p0 = [-2.1, 0.4, -12.8]
    p = p0
    cnt=0
    for i in arange(-12,9,3.0):
        p[0] = i
        ax1.plot(xx, cnt+rz.dtdfunc(xx, *p, norm=False), 'k-',alpha = 0.3)
        cnt +=0.1 

    p = p0
    cnt = 0
    for i in arange(0.1,12.,0.5):
        p[1] = i
        ax2.plot(xx, cnt+rz.dtdfunc(xx, *p), 'k-',alpha = 0.3)
        cnt +=0.1 

    p = p0
    cnt = 0
    for i in arange(-12,9,3.0):
        p[2] = i
        ax3.plot(xx, cnt+rz.dtdfunc(xx, *p), 'k-',alpha = 0.3)
        cnt+=0.1


    show()
    
             
    
