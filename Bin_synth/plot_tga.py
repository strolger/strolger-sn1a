#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit

if __name__=='__main__':



    ax=subplot(111)
    ax2=ax.twinx()
    data=loadtxt('toonan_ga.txt')


    p1 = [-1257.93, 59.32, 248.88] ## from MCMC
    time = arange(0.05,13.5,0.1)
    scale = 12.921*(0.0210)*0.062*(0.7)**2
    dtd = scale*rz.dtdfunc(time,*p1)

    def fn(x,*p):
        aa=p
        return(rz.dtdfunc(x,*p1)*scale*aa)
    
    p0 = [1.0]
    popt,pcov = curve_fit(fn, data[:,0],data[:,2],p0=p0)
    print(popt)
    ax.plot(time,dtd,'b-', lw=2, label= 'Exponential model')#label='Norm = %1.1f' %(simps(dtd,x=time)))
    ax2.plot(data[:,0], data[:,2]/(popt*0.5),'k-',label=r'$\alpha\alpha$ model, Toonan et al. 2012')

    pwrl = (-1.39, 1.0)
    dtd = rz.powerdtd(time, *pwrl, normed=False)*5.4e-3

    def fn2(x,*p):
        aa = p
        return(rz.powerdtd(x, *pwrl, normed=False)*5.4e-3*aa)
     
    p0 = [1.0]
    popt,pcov = curve_fit(fn2, data[:,0],data[:,1],p0=p0)
    print(popt)
   
    ax.plot(time, dtd, 'r-', label=r'Power-law model; $\beta=-1$')
    ax2.plot(data[:,0], data[:,1]/popt,'k--',label=r'$\gamma\alpha$ model')
    
    ax.set_yscale('log')
    ax2.set_yscale('log')
    ax.set_xlabel('Delay time (Gyr)')
    ax.set_ylabel(r'SN Ia yr$^{-1}$ ($10^{10}$ M$_{\odot}$)$^{-1}$')
    ax2.set_ylabel(r'Merger rate  yr$^{-1}$ ($10^{10}$ M$_{\odot}$)$^{-1}$')

    ax.set_ylim(2e-5, 0.1)
    ax2.set_ylim(2e-5, 0.1)

    ax.legend(frameon=False)
    ax2.legend(loc=3,frameon=False)
    savefig('figure_toonan.png')
    
