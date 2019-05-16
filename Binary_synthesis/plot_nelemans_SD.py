#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit


if __name__=='__main__':


    data = loadtxt('nelemans_table.txt')
    #data = data[:,:-3]
    ax=subplot(111)
    tt = arange(0,10000,100)
    p0 = [
        [10, 600, 220],
        [110, 1000, 2],
        [350, 1200, 20],
        [6000, 6000, -2],
        ]
    scales = [
        6.,
        2.,
        6.,
        8e-2,
        ]
        
    ## for ii in range(shape(data)[1]-2):
    ##     ax.errorbar(data[:,0], data[:,ii+2]*1e-13, xerr=data[:,1], fmt='o')
    ##     ax.plot(tt, rz.dtdfunc(tt, *p0[ii])*scales[ii])

    ii=0
    ax.errorbar(data[:,0], data[:,2]*1e-3, xerr=data[:,1], fmt='o',color=u.my_color(0), label='Mennekens et al. 2010')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(0))
    ii=1
    ax.errorbar(data[:,0], data[:,3]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(2), label='Yungelson 2010')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(2))
    ii=2
    ax.errorbar(data[:,0], data[:,5]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(5), label='Wang et al. 2010')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(5))
    ii=3
    ax.errorbar(data[:,0], data[:,6]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(6), label='Ruiter et al. 2009')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(6))

    

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(4e-7,2e-2)
    ax.legend(loc=3)

    ax.set_xlabel('Delay time (Myr)')
    ax.set_ylabel(r'SN Ia yr$^{-1}$ ($10^{10}$ M$_{\odot}$)$^{-1}$')
    savefig('figure_sd_fits.png')
    
    
    
