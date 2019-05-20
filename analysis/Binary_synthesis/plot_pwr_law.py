#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit


if __name__=='__main__':


    fig= plt.figure()
    grid = plt.GridSpec(1,3,wspace=0.2)

    data = loadtxt('nelemans_table.txt')

    ax = fig.add_subplot(grid[0,:-1])
    ax2 = fig.add_subplot(grid[0,2])
    tt = arange(0,10000,100)
        
    ax.errorbar(data[:,0], data[:,2]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(0), label='Mennekens et al. 2010')
    ax.errorbar(data[:,0], data[:,3]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(2), label='Yungelson 2010')
    ax.errorbar(data[:,0], data[:,5]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(5), label='Wang et al. 2010')
    ax.errorbar(data[:,0], data[:,6]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(6), label='Ruiter et al. 2009')

    scale = (0.0210)*0.062*(0.7)**2
    p1 = (-1258, 59, 248)
    ax.plot(tt,5.28*scale*rz.dtdfunc(tt/1e3, *p1),'k-') 

    files = glob.glob('../DTD_fits/mc_sfd_*.pkl')
    tot = 0
    for file in files:
        samples = pickle.load(open(file,'rb'))
        samples = samples[:,50:,:].reshape((-1,5))
        print('adding %d samples from %s... '%(len(samples), file))
        try:
            temp = concatenate((temp, samples), axis=0)
        except:
            temp = samples

    print('%d total samples' %len(temp))
    samples = temp
    ndraws = 100 #21 is actually pretty good
    draws  = samples[np.random.randint(0, len(samples), ndraws),:]
    ysample = asarray([rz.dtdfunc(tt/1e3, *pi[1:-1])*exp(pi[0])*(0.0210)*(0.7)**2*5.28 for pi in draws])


    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(tt, upper, lower, color='green', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(tt, upper, lower, color='green', alpha=0.2)



    ## ax2 =subplot(122)
    data = loadtxt('nelemans_table2.txt')
    ax2.errorbar(data[:,0], data[:,2]*1e-3, xerr=data[:,1], fmt='o',color=u.my_color(0), label='Yungelson 2010')

    pa = (-1.0, 1.0)
    yy = rz.powerdtd(tt, *pa)
    ax2.plot(tt, yy*2.0, color=u.my_color(2), label=r'$t^{-1}$')
    ## p1=(-662, 2200, 1101)
    ## ax2.plot(tt, 5*rz.dtdfunc(tt,*p1), '--', color=u.my_color(5), label='best model fit')



    ax2.plot(tt,5.28*scale*rz.dtdfunc(tt/1e3, *p1),'k-') 
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax2.fill_between(tt, upper, lower, color='green', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax2.fill_between(tt, upper, lower, color='green', alpha=0.2)
    

    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax2.set_xlim(100,10000)
    ## ax2.set_ylim(1e-5, 0.005)
    ax2.set_ylim(4e-7, 2e-2)
    ax2.set_xlabel('Delay time (Myr)')


    ax2.set_yticks([])
    ax2.legend(loc=3, fontsize=9,frameon=False)
    ax2.set_title('DD Models')

    

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(4e-7,2e-2)
    ax.legend(loc=3, frameon=True)

    ax.set_xlabel('Delay time (Myr)')
    ax.set_ylabel(r'SN Ia yr$^{-1}$ ($10^{10}$ M$_{\odot}$)$^{-1}$')
    ax.set_title('SD Models')
    savefig('figure_pwr_law.png')
    
    
    
