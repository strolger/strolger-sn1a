#!/usr/bin/env python
'''
this wont work as a straight least squares fit
But does favor more in the delayed bin than in the
prompt bin.
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


def step_dtd(x,*p, norm=False):
    A,B, C=p
    ret = zeros((len(x)),)
    ret[where((0.0< x) & (x<= 1.0))]=A
    ret[where((1.0< x) & (x<= 5.0))]=B
    ret[where((5.0< x) & (x<= 13.6))]=C
    return(ret)

def sdtdfit(time,*p):
    ff, aa, bb, cc = p
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale *0.7**2.*1e4
    sfh = rz.sfr_behroozi_12(time)
    dt = sum(diff(time))/(len(time)-1)
    p1 = (aa, bb, cc)
    res = step_dtd(time,*p1)
    tt = arange(0,13.6,0.01)
    norm_res = sum(step_dtd(tt,*p1))*0.01
    res = res/norm_res
    tmp = convolve(sfh, res, 'full')
    return(ff*tmp[:len(time)]*dt*scale)

def fit_one(rates, *p):
    data = deepcopy(rates)
    tt = 13.6-array([ct.cosmotime(x) for x in data[:,0]]) ## now also in forward time
    data[:,0] = tt
    data = data[argsort(data[:,0])]
    p0=(0.05,)+p
    popt, pcov = curve_fit(sdtdfit, data[:,0], data[:,1], p0=p0,
                           sigma=data[:,2],
                           bounds=(0,[1.,100.,100.,100.]))#, method='dogbox')
    return(popt,pcov)


def plot_one(rates,plotname,*p, frac=0.05, age=13.6):
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4 ## factors of h...
    dt = 0.05
    tt = arange(0.1,age,dt) ## forward time
    lbt = age - tt
    zz = [ct.cosmoz(x) for x in lbt]

    sfh = rz.sfr_behroozi_12(tt)
    dtd = step_dtd(tt, *p)
    norm_dtd = sum(step_dtd(tt,*p))*dt
    dtd = dtd/norm_dtd
    
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
    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'SN Ia Rate')
    ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ax2.set_ylabel(r'SF Rate')
    
    ax3 = axes([0.65, 0.6, 0.23, 0.2])
    ax3.plot(tt,dtd,'b-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
    ax3.plot(tt,rz.powerdtd(tt, *pwrl), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))
    pn = (-20,8.8, 11)
    ax3.plot(tt,rz.dtdfunc(tt,*pn),'b--',label='Best model')
    ax.set_title(r' $k=%.4f$ M$_{\odot}^{-1}$, $f=%2.1f\%%$' %(scale_k,frac*100))
 

    ax3.set_ylabel('$\Phi$')
    ax3.set_xlabel('Delay Time (Gyr)')
    ax3.set_xlim(0,12)
    ax3.set_ylim(0,1.3)
    ax3.legend(frameon=False)
    savefig(plotname)
    return()
    

if __name__=='__main__':

    p0 = (10., 1.,1.)
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())
    popt,pcov=fit_one(rates,*p0)
    print(popt, pcov)
    print(sqrt(diag(pcov))[0],sqrt(diag(pcov))[1])#,sqrt(diag(pcov))[2])
    plot_one(rates,'figure_sfd_optimized_bins.png',*popt[1:], frac=0.05)
    ## plot_one(rates,'figure_sfd_optimized_bins.png',*p0, frac=0.07)

    
    
