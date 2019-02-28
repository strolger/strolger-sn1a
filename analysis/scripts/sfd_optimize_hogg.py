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
import warnings
warnings.simplefilter("ignore",RuntimeWarning)



def dtdfit(time,*p):
    ff, aa, bb, cc = p
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale *0.7**2.*1e4
    par_model = [0.0134, 2.55, 3.3, 6.1]
    sfh = rz.csfh_time(time, *par_model)
    dt = sum(diff(time))/(len(time)-1)
    p1 = (aa, bb, cc)
    res = rz.dtdfunc(time, *p1,norm=True)
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
    zz = [ct.cosmoz(x, ho=70) for x in lbt]

    par_model = [0.0134, 2.55, 3.3, 6.1]

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
    ax2.plot(zz, sfh, 'r-', label='CSFH')
    ax.plot(zz, rate_fn, 'k-', label='Fit')
    ax.plot(zz, jdud, 'k:', label=r'$t^{%.1f}$'%(pwrl[0]))
    
    ax.errorbar(rates[:,0], rates[:,1], yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6', alpha=0.4)
    ax.errorbar(brates[:,0], brates[:,1], yerr=[brates[:,3],brates[:,2]],
                xerr=[brates[:,5],brates[:,4]], fmt='o', color='0.0',zorder=10)

    pwrl = (-1.0,1.0)
    ax.set_xlim(0,2.6)
    ax.set_ylim(0,1.8)
    ax2.set_ylim(0,0.16)
    
    ax3 = axes([0.65, 0.62, 0.23, 0.2])
    ax3.plot(tt,log10(dtd),'b-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
    ax3.plot(tt,log10(rz.powerdtd(tt, *pwrl)), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))
    ax3.set_ylabel('$\log(\Phi)$')
    ax3.set_xlabel('Delay Time (Gyr)')
    ax3.set_xlim(0,12)
    ax3.set_ylim(-3,0.5)
    ##ax3.set_ylim(0,1)
    ## ax3.legend(loc=1,ncol=1,frameon=False, fontsize=8)
    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ax2.set_ylabel(r'M$_{\odot}$ per year per Mpc$^3$')
    ## ax.set_title(r' $k=%.4f$ M$_{\odot}^{-1}$, $f=%2.1f\%%$' %(scale_k,frac*100))
    ax.set_title(r'$k=%.4f\,M_{\odot}^{-1},\,\,\varepsilon=%2.1f\%%$' %(scale_k,frac*100))
    ax.legend(loc=2, frameon=False)
    ax2.legend(loc=4, frameon=False)
    
    savefig(plotname)
    return()




if __name__=='__main__':

    calculate = False
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,verbose=False,splits=arange(0,1.167,0.167).tolist())
    data = deepcopy(rates) #this seems to work better than then actual values
    #data = deepcopy(rates[rates[:,0].argsort()])
    if calculate:
        tt = 13.6-array([ct.cosmotime(x, ho=70) for x in data[:,0]])
        data[:,0] = tt
        data = data[argsort(data[:,0])]

        import scipy.optimize as op
        nll = lambda *args: -lnlike(*args)
        p0 = (0.05, 3.5, 2.5, 2.5, 0.01)
        res = op.minimize(nll, [p0[0],p0[1],p0[2],p0[3],log(p0[4])],
                          bounds=[(0.,1.0),(-2000.,2000.),(0.001,100.),(-500.0,500.0),(-4.,0.)],
                          args=(data[:,0], data[:,1], data[:,3]))

        p2=res['x'][1:4]
        frac = res['x'][0]
        print(p2,frac)

    else:
        p2 = [-1518.39108091,    51.06020462,    49.98728772] #0.0976723208669[-593.54527392,   41.46938593,  173.52397092]
        frac = 0.062
    plot_one(rates,'figure_sfd_optimized.png',*p2, frac=frac)
    
