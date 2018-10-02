#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import cosmotools as ct
from scipy.integrate import simps,quad
from strolger_util import imf
from strolger_util import forkmap


parallel=False
Nproc=10

def run_model(tt, *p):
    ddt = average(diff(tt))
    csfh_model = [0.013, 2.6, 3.2, 6.1]
    sfh = rz.csfh_time(tt, *csfh_model)
    dtd = rz.dtdfunc(tt, *p)
    tmp = convolve(sfh, dtd, 'full')*ddt#*frac*scale ## now convolved result in forward time
    rate_fn=tmp[:len(dtd)]
    print(p,sum(rate_fn))
    if sum(rate_fn)==0.:
        pdb.set_trace()
    return(rate_fn)
    

if __name__=='__main__':

    dt = 0.05
    age=13.6
    tt = arange(dt,age,dt)
    p_val = [-1031, 65, 190]
    p_err = [[1033., -687.],
             [15., -20.],
             [210., -210.]
             ]


    ax = subplot(121)
    ## ax3 = ax.twinx()
    
    ## ax2 = axes([0.66, 0.65, 0.23, 0.2])
    ax2 = subplot(122)
    ax2.plot(tt, rz.dtdfunc(tt, *p_val), 'b-')

    p_err=array(p_err)
    p_err=sqrt(p_err[:,0]**2.+p_err[:,1]**2.)


    tta = linspace(min(tt), max(tt), 100)
    cov = diag(array(p_err)**2.)
    ps = np.random.multivariate_normal(p_val, cov, 100)
    ## ps = np.random.multivariate_normal(p_val, cov, 10000)
    ysample = asarray([rz.dtdfunc(tta, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)

    ax2.fill_between(tta, upper, lower, color='green', alpha=0.2)
    ## ax.plot(tta, upper, 'b--')
    ## ax.plot(tta, lower, 'b:')

    pwrl=(-1.,1.)
    ax2.plot(tt,rz.powerdtd(tt, *pwrl), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))

    
    ax2.set_xlim(0,12)
    ax2.set_ylim(0,1.1)

    frac = 0.062
    ## p =  (-1.12659790e+03,   5.69300820e+01,   1.20864134e+02)
    ## age = 13.6
    
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.167,0.167).tolist())
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4## factors of h...
    lbt = age - tta
    zz = [ct.cosmoz(x, ho=70) for x in lbt]


    ax.errorbar(rates[:,0], rates[:,1], yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6', alpha=0.4)
    ax.errorbar(brates[:,0], brates[:,1], yerr=[brates[:,3],brates[:,2]],
                xerr=[brates[:,5],brates[:,4]], fmt='o', color='0.0',zorder=10)

    ###[NOTE: need a way to run through convolutions as a function...]

    if not parallel:
        ysample = asarray([run_model(tta, *pi) for pi in ps])*frac*scale
    else:
        @forkmap.parallelizable(Nproc)
        def run_ysample(*p):
            print(p)
            yi = run_model(tta, *p)*frac*scale
            print (p, yi)
            return(p, yi)
        ysample = forkmap.map(run_ysample, *ps)
        pdb.set_trace()

    ## a = ones((len(ps),1),)*frac
    ## ps = hstack((a,ps))

    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)

    ax.plot(zz, run_model(tta, *p_val)*frac*scale, 'b-')
    ## ax.plot(zz, upper, 'b--')
    ## ax.plot(zz, lower, 'b:')
    ax.fill_between(zz, upper, lower, color='green', alpha=0.2)
    ax.set_xlim(0,2.5)

    savefig('temp2.png')
    sys.exit()




    ## par_model = [0.013, 2.6, 3.2, 6.1]
    ## sfh = rz.csfh_time(tt, *par_model)
    ## dtd = rz.dtdfunc(tt, *p)
    
    ## tmp = convolve(sfh, dtd, 'full')*dt*frac*scale ## now convolved result in forward time
    ## rate_fn=tmp[:len(dtd)]
    ## clf()
    ## ax = subplot(111)
    ## ax2 = ax.twinx()
    ## ax2.plot(zz, sfh, 'r-')
    ## ax.plot(zz, rate_fn, 'k-')
    ## ax.axhline(0., color='k')
    ## ax.errorbar(rates[:,0], rates[:,1], yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6', alpha=0.4)
    ## ax.errorbar(brates[:,0], brates[:,1], yerr=[brates[:,3],brates[:,2]],
    ##             xerr=[brates[:,5],brates[:,4]], fmt='o', color='0.0',zorder=10)
   

    
    ## #junk, yn = u.recast(rates[:,0], 0., zz, rate_fn)
    ## #ax.errorbar(rates[:,0], rates[:,1] - yn, yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6')
    

    ## pwrl = (-1.0,1.0)
    ## ax.set_xlim(0,2.5)
    ## ax.set_ylim(0,1.8)

    ## ax2.set_ylim(0,0.16)
    
    ## ax3 = axes([0.66, 0.65, 0.23, 0.2])
    ## ax3.plot(tt,dtd,'b-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
    ## ax3.plot(tt,rz.powerdtd(tt, *pwrl), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))
    ## ax3.set_ylabel('$\Phi$')
    ## ax3.set_xlabel('Delay Time (Gyr)')
    ## ax3.set_xlim(0,12)
    ## ax3.set_ylim(0,1.3)
    ## ax3.legend(frameon=False)

    ## ax.set_xlabel('Redshift')
    ## ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ## ax2.set_ylabel(r'M$_{\odot}$ per year per Mpc$^3$')
    ## ## ax.set_title(r' $k=%.4f$ M$_{\odot}^{-1}$, $f=%2.1f\%%$' %(scale_k,frac*100))
    ## ax.set_title(r'$k=%.4f\,M_{\odot}^{-1},\,\,f=%2.1f\%%$' %(scale_k,frac*100))
    ## savefig('temp.png')
