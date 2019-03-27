#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit

if __name__=='__main__':


    ## mock data
    tt = arange(0.04,13.6,0.1)
    p_pwr = (-1., 1.)
    yy = rz.powerdtd(tt, *p_pwr)#, normed=False)

    ax=subplot(111)
    ax.plot(tt, yy, 'r--')
    #ax.set_xlim(0.04,14)


    p0=(2,2,1)
    popt,pcov=curve_fit(rz.dtdfunc, tt, yy, p0=p0)
    ax.plot(tt, rz.dtdfunc(tt,*popt), '-',color='orange')
    print(popt)
    p1 = (0.15, 2, -0.7)
    ax.plot(tt, rz.dtdfunc(tt,*p1), 'r-')
    
    savefig('junk.png')
    
              
