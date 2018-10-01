#!/usr/bin/env python
'''
Gives an optimal fit to binned rate data
using curve fit
L. Strolger
2018
'''
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import imf
from strolger_util import cosmotools as ct
from scipy.integrate import simps,quad
from scipy.optimize import curve_fit
from copy import copy, deepcopy
import warnings
warnings.simplefilter("ignore",RuntimeWarning)

def dtdfit(time,*p):
    ff, aa, bb, cc = p
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale *0.7**2.*1e4
    par_model = [0.013, 2.6, 3.2, 6.1]
    sfh = rz.csfh_time(time, *par_model)
    dt = sum(diff(time))/(len(time)-1)
    p1 = (aa, bb, cc)
    res = rz.dtdfunc(time, *p1,norm=True)
    tmp = convolve(sfh, res, 'full')
    return(ff*tmp[:len(time)]*dt*scale)

def fit_one(rates, *p):
    data = deepcopy(rates)
    tt = 13.6-array([ct.cosmotime(x, ho=70) for x in data[:,0]])
    data[:,0] = tt
    data = data[argsort(data[:,0])]    
    p0 = (0.05,)+p
    print (p0)
    popt, pcov = curve_fit(dtdfit, data[:,0], data[:,1], p0=p0,
                           sigma=data[:,3],
                           ## bounds=([0.,-1.0e5,0.001,-1000.],[1.,1e5,1000.,1000.]),
                           maxfev=10000)
    return(popt,pcov)

def plot_one(rates,plotname,*p, frac=0.05, age=13.6):
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4 ## factors of h...
    dt = 0.05
    tt = arange(0.0,age,dt)
    lbt = age - tt
    zz = [ct.cosmoz(x, ho=70) for x in lbt]

    par_model = [0.013, 2.6, 3.2, 6.1]
    sfh = rz.csfh_time(tt, *par_model)
    dtd = rz.dtdfunc(tt, *p)
    
    tmp = convolve(sfh, dtd, 'full')*dt*frac*scale ## now convolved result in forward time
    rate_fn=tmp[:len(dtd)]
    clf()
    ax = subplot(111)
    ax2 = ax.twinx()
    ax2.plot(zz, sfh, 'r-')
    ax.plot(zz, rate_fn, 'k-')
    
    ax.errorbar(rates[:,0], rates[:,1], yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6')
    ax.errorbar(brates[:,0], brates[:,1], yerr=[brates[:,3],brates[:,2]],
                xerr=[brates[:,5],brates[:,4]], fmt='o', color='0.0',zorder=10)

    pwrl = (-1.0,1.0)
    ax.set_xlim(0,2.5)
    ax.set_ylim(0,1.8)
    ax2.set_ylim(0,0.16)
    ax3 = axes([0.65, 0.6, 0.23, 0.2])
    ax3.plot(tt,dtd,'b-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
    ax3.plot(tt,rz.powerdtd(tt, *pwrl), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))
    ax3.set_ylabel('$\Phi$')
    ax3.set_xlabel('Delay Time (Gyr)')
    ax3.set_xlim(0,12)
    ax3.set_ylim(0,1.3)
    ax3.legend(frameon=False)

    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ax2.set_ylabel(r'M$_{\odot}$ per year per Mpc$^3$')
    ax.set_title(r'$k=%.4f\,M_{\odot}^{-1},\,\,f=%2.1f\%%$' %(scale_k,frac*100))

    savefig(plotname)
    return()
    

if __name__=='__main__':

    p0 = (0.5, 0.5, 2.2)
    ## p0 = (-1.14697256e+03, 6.02664905e+01,  1.15263257e+02)
    
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())#, verbose=True)
    popt,pcov=fit_one(brates,*p0)
    print(popt,pcov)
    print(sqrt(diag(pcov))[0],sqrt(diag(pcov))[1],sqrt(diag(pcov))[2],sqrt(diag(pcov))[3])
    plot_one(rates,'figure_sfd_optimized_curvefit.png',*popt[1:],frac=popt[0])
    
    
    
