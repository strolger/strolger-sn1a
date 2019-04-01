#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit

if __name__=='__main__':


    ## mock data
    tt = arange(0.04,13.6,0.1)
    p_pwr = (-1, 1.)
    yy = rz.powerdtd(tt, *p_pwr)#, normed=False)



    zz = yy*0.0
    id1 = where(tt < 0.5)
    zz[id1]=sum(yy[id1])
    id2 = where(tt >= 0.5)
    zz[id2]=sum(yy[id2])

    
    ax=subplot(111)
    ax.plot(tt, yy, 'r--')
    ax.plot(tt,zz, 'b--')
    p0 = (1,1,1)
    popt,pcov=curve_fit(rz.dtdfunc, tt, yy, p0=p0)
    ax.plot(tt, rz.dtdfunc(tt,*popt), '-',color='orange')
    print(popt,pcov)
    p1=(-1258, 59, 248)

    yy = rz.dtdfunc(tt,*p1)

    zz = yy*0.0
    id1 = where(tt < 0.5)
    zz[id1]=sum(yy[id1])
    id2 = where(tt >= 0.5)
    zz[id2]=sum(yy[id2])

    ax.plot(tt,yy , 'r-')
    ax.plot(tt,zz , 'b-')


    
    #ax.set_xlim(0.04,14)
    savefig('junk.png')
    
              
