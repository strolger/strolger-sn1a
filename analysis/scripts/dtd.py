#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz






if __name__=='__main__':

    xx = arange(0,10,0.1)

    ax1 = subplot(131)
    ax2 = subplot(132)
    ax3 = subplot(133)

    p = (5,5,5)

    for i in arange(-12,9,3.0):
        p = (i, 5, 5)
        ax1.plot(xx, rz.dtdfunc(xx, *p), 'k-',alpha = 0.3)

    for i in arange(1.0,12.,1.0):
        p = (5, i, 5)
        ax2.plot(xx, rz.dtdfunc(xx, *p), 'k-',alpha = 0.3)

    for i in arange(-12,9,3.0):
        p = (5, 5, i)
        ax3.plot(xx, rz.dtdfunc(xx, *p), 'k-',alpha = 0.3)

    show()
    
             
    
