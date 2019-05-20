#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import cosmotools as ct
from scipy.integrate import simps,quad
from strolger_util import imf

csfh_model = [0.0134, 2.55, 3.3, 6.1]


def run_model(tt, *p):
    ddt = average(diff(tt))
    sfh = rz.csfh_time(tt, *csfh_model)
    dtd = rz.dtdfunc(tt, *p[1:])
    tmp = convolve(sfh, dtd, 'full')*ddt*np.exp(p[0])#*frac*scale ## now convolved result in forward time
    rate_fn=tmp[:len(dtd)]
    return(rate_fn)

    

if __name__=='__main__':

    dt = 0.05
    age=13.6
    tt = arange(dt,age,dt)
    ## p_val = [-2.810, -1518.39108091,    51.06020462,    49.98728772] ## from optimized fit
    p_val = [-2.814, -1257.93, 59.32, 248.88] ## from MCMC

    ax = subplot(111)
    ax2 = axes([0.55, 0.6, 0.33, 0.25])
    ## ax2 = axes([0.65, 0.62, 0.23, 0.2])
    ax2.plot(tt, log10(rz.dtdfunc(tt, *p_val[1:])), 'b-')
    
    files = glob.glob('mc_sfd_*.pkl')
    tot = 0
    for file in files:
        samples = pickle.load(open(file,'rb'))
        samples = samples[:,50:,:].reshape((-1,5))
        print('adding %d samples from %s... '%(len(samples), file))
        try:
            temp = concatenate((temp, samples), axis=0)
        except:
            temp = samples
        
    #tta = linspace(min(tt), max(tt), 100)
    tta = tt
    print('%d total samples' %len(temp))
    samples = temp
    ndraws = 100 #21 is actually pretty good
    draws  = samples[np.random.randint(0, len(samples), ndraws),:]
    ysample = asarray([rz.dtdfunc(tta, *pi[1:-1]) for pi in draws])


    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax2.fill_between(tta, log10(upper), log10(lower), color='green', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax2.fill_between(tta, log10(upper), log10(lower), color='green', alpha=0.2)

    pwrl=(-1.,1.)
    ax2.plot(tt,log10(rz.powerdtd(tt, *pwrl)), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))

    
    frac = 0.065
    
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,splits=arange(0,1.125,0.125).tolist())
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4## factors of h...
    lbt = age - tta
    zz = [ct.cosmoz(x, ho=70) for x in lbt]


    ax.errorbar(rates[:,0], rates[:,1], yerr=[rates[:,3],rates[:,2]], fmt='o', color='0.6', alpha=0.4)
    ax.errorbar(brates[:,0], brates[:,1], yerr=[brates[:,3],brates[:,2]],
                xerr=[brates[:,5],brates[:,4]], fmt='o', color='0.0',zorder=10)

    ysample = asarray([run_model(tta, *pi[:-1]) for pi in draws])*scale
        

    ax.plot(zz, run_model(tta, *p_val)*scale, 'b-', label='Fit')
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(zz, upper, lower, color='green', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(zz, upper, lower, color='green', alpha=0.2)


    pwrl = (-1.0,1.0)
    dud = rz.powerdtd(tta, *pwrl)
    jdud = convolve(rz.csfh_time(tta, *csfh_model), dud, 'full')*average(diff(tta))*frac*scale
    jdud = jdud[:len(dud)]
    ax.plot(zz, jdud, 'b:', label = r'$t^{%.1f}$'%(pwrl[0]))
    ## pdb.set_trace()

    ax.set_xlim(0,2.5)
    ax.set_ylim(0,1.8)
    ax2.set_xlim(0,12)
    ax2.set_ylim(-3,0.5)

    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ax2.set_ylabel('$\log(\Phi)$')
    ax2.set_xlabel('Delay Time (Gyr)')
    ax.set_title(r'$k=%.4f\,M_{\odot}^{-1},\,\,\varepsilon=%2.1f\%%$' %(scale_k, np.exp(p_val[0])*100))
    ax.legend(loc=2, frameon=False)
    savefig('figure_fit_demo_werr.png', transparent=True)



