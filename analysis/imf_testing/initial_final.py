#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u


def m_f(x):
    yy = 0.187*x+0.184
    yy[where(x < 2.85)]=0.08*x[where(x<2.85)]+0.489
    yy[where(x > 3.60)]=0.107*x[where(x>3.60)]+0.471
    return(yy)
def m_f1(x):
    yy = 0.126*x-0.015
    yy[where(x < 2.85)]=0.064*x[where(x<2.85)]+0.459
    yy[where(x > 3.60)]=0.123*x[where(x>3.60)]+0.400
    return(yy)
def m_f2(x):
    yy = 0.248*x+0.383
    yy[where(x < 2.85)]=0.096*x[where(x<2.85)]+0.519
    yy[where(x > 3.60)]=0.123*x[where(x>3.60)]+0.548
    return(yy)


mi = arange(0,8,0.1)


ax = subplot(111)
ax.plot(mi, m_f(mi), 'k-')
ax.plot(mi, m_f1(mi), 'r--')
ax.plot(mi, m_f2(mi), 'r--')

inn = array([0.6,0.65,0.7])
junk, out= u.recast(inn, 0., m_f1(mi), mi)
print(out)

show()
