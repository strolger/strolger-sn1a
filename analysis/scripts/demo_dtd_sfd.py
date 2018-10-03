#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import imf
from strolger_util import cosmotools as ct
from scipy.integrate import simps,quad
from copy import copy, deepcopy



def plot_one(rates,plotname,*p, frac=0.05, age=13.6):
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4 ## factors of h...
    dt = 0.05
    tt = arange(dt,age,dt) ## forward time
    lbt = age - tt
    zz = [ct.cosmoz(x) for x in lbt]

    par_model = [0.013, 2.6, 3.2, 6.1]
    sfh = rz.csfh_time(tt, *par_model)
    dtd = rz.dtdfunc(tt, *p)
    pwrl = (-1.0,1.0)
    dud = rz.powerdtd(tt, *pwrl)
    
    tmp = convolve(sfh, dtd, 'full')*dt*frac*scale ## now convolved result in forward time
    rate_fn=tmp[:len(dtd)]

    jdud = convolve(sfh, dud, 'full')*dt*scale*0.065
    jdud = jdud[:len(dtd)]
    clf()
    ax = subplot(111)
    ax2 = ax.twinx()
    ax2.plot(zz, sfh, 'r-')
    ax.plot(zz, rate_fn, 'k-')
    ax.plot(zz, jdud, 'k:')
    
    ax.errorbar(rates[:,0], rates[:,1], yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6')
    ax.errorbar(brates[:,0], brates[:,1], yerr=[brates[:,3],brates[:,2]],
                xerr=[brates[:,5],brates[:,4]], fmt='o', color='0.0',zorder=10)

    
    ax.set_xlim(0,2.5)
    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'SN Ia Rate')
    ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ax2.set_ylabel(r'SF Rate')
    
    ax3 = axes([0.65, 0.6, 0.23, 0.2])
    ax3.plot(tt,dtd,'b-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
    ax3.plot(tt,dud, 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))
 

    ax3.set_ylabel('$\Phi$')
    ax3.set_xlabel('Delay Time (Gyr)')
    ax3.set_xlim(0,12)
    ax3.set_ylim(0,1.3)
    ax3.legend(frameon=False)
    savefig(plotname)
    return()

if __name__=='__main__':
    
    ## p0 = (0.05, -5.1, 6.5, 2.1)
    ## p0 = (0.1, 3.3, 7.1, 0.7)
    ## p0 = (0.05, -2.99, 12.61, -7.6)
    p0 = (0.05, 3.5, 0.5, 2.2)
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    plot_one(rates,'figure_fit_demo.png',*p0[1:], frac=p0[0])

    
