#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit

if __name__=='__main__':

    tt = arange(-4000, 10000, 10)


    ## dd = loadtxt('maoz_data.txt')
    ## ## dd = loadtxt('maoz_table2.txt')
    ax = subplot(111)


    ## ax.errorbar(dd[:,0],dd[:,3],
    ##             xerr=[abs(dd[:,1]), dd[:,2]],
    ##             yerr=[abs(dd[:,4]), dd[:,5]],
    ##             fmt='o', color='k')

    data = loadtxt('nelemans_table2.txt')
    ax.errorbar(data[:,0], data[:,2]*1e-3, xerr=data[:,1], fmt='o',color=u.my_color(0), label='Mennekens et al. 2010')

    pa = (-1.0, 1.0)
    yy = rz.powerdtd(tt, *pa)
    ax.plot(tt, yy*0.5)
    p1=(-662, 2200, 1101)
    ax.plot(tt, 5*rz.dtdfunc(tt,*p1), '--')


    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(100,10000)
    ax.set_ylim(1e-5, 0.005)
    savefig('figure_pwr.png')
    
    
               
