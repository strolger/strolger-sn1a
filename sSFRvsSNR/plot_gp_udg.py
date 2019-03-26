#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u

def func1(x, *p):
    d,tau = p
    return(x**d*exp(-x/tau))

def func2(x,*p):
    B,C,tau=p
    return(((x/tau)**B + (x/tau)**(-C))**(-1))




if __name__=='__main__':



    
    xx = arange(0.05, 13.6, 0.005)
    p0a = (7.6, 0.5)
    p0b = (18.6, 0.8, 6.7)

    ax = subplot(111)
    ax.plot(xx,func1(xx, *p0a),'r-')
    ## ax.plot(xx,14*func2(xx, *p0b),'r--')
    ## ax.plot(xx, func1(xx, *p0a)+14*func2(xx,*p0b), 'b-')

    show()
