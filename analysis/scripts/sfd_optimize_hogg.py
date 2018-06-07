#!/usr/bin/env python
'''
Gives an optimal fit to binned rate data
using Hoggs generative model
https://arxiv.org/pdf/1008.4686.pdf
L. Strolger
2018
'''
import os,sys,pdb,scipy,glob,pickle,datetime
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import imf
from strolger_util import cosmotools as ct
from copy import copy, deepcopy
from scipy.integrate import simps,quad
import emcee

def dtdfit(time,*p):
    ff, m, w, k = p
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale *0.7**2.*1e4
    sfh = rz.sfr_behroozi_12(time)
    dt = sum(diff(time))/(len(time)-1)
    p1 = (m, w, k)
    res = rz.dtdfunc(time, *p1)
    tmp = convolve(sfh, res, 'full')
    return(ff*tmp[:len(time)]*dt*scale)


def lnlike(p, x, y, yerr):
    ff, m, w, k, lnf = p
    p1 = (ff, m, w, k)
    model = dtdfit(x,*p1)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))


def plot_one(rates,plotname,*p, frac=0.05, age=13.6):
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4## factors of h...
    dt = 0.05
    tt = arange(dt,age,dt)
    lbt = age - tt
    zz = [ct.cosmoz(x) for x in lbt]

    sfh = rz.sfr_behroozi_12(tt)
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
    ax.set_title(r' $k=%.4f$ M$_{\odot}^{-1}$, $f=%2.1f\%%$' %(scale_k,frac*100))
    
    savefig(plotname)
    return()




if __name__=='__main__':

    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,verbose=True,splits=arange(0,1.167,0.167).tolist())
    data = deepcopy(brates) #this seems to work better than then actual values
    tt = 13.6-array([ct.cosmotime(x) for x in data[:,0]])
    data[:,0] = tt
    data = data[argsort(data[:,0])]
    
    import scipy.optimize as op
    nll = lambda *args: -lnlike(*args)
    p0 = (0.05, 3.5, 2.5, 2.5, 0.5)
    res = op.minimize(nll, [p0[0],p0[1],p0[2],p0[3],log(p0[4])],
                      bounds=[(0.,1.0),(-2000.,2000.),(0.001,100.),(-500.0,500.0),(-4.,0.)],
                      args=(data[:,0], data[:,1], data[:,3]))
    
    p2=res['x'][1:4]
    
    print(res['x'])
    plot_one(rates,'figure_sfd_optimized.png',*p2, frac=res['x'][0])
    pdb.set_trace()
