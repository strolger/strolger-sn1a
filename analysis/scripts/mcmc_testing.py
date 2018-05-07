#!/usr/bin/env python
'''
planying, and learning, emcee using examples
on
http://dfm.io/emcee/current/user/line/#maximum-likelihood-estimation
'''

import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
import emcee


def linef(x,*p):
    m,B = p
    return(m*x+B)



def synthetic_data (*p):
    m_true, b_true, f_true=p
    N = 50
    x = sort(10*np.random.rand(N))
    yerr = 0.1+0.5*np.random.rand(N)
    p1 = (m_true, b_true)
    y = linef(x,*p1)
    y+=abs(f_true*y) * np.random.randn(N)
    y+=yerr * np.random.randn(N)
    return(x,y,yerr)



def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = linef(x,*(m,b)) 
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))


def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf


def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)



if __name__=='__main__':


    showburn=False
    showcorner=False
    
    trueval = (-0.95, 4.29, 0.53)
    x, y, yerr = synthetic_data(*trueval)
    import scipy.optimize as op
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [trueval[0], trueval[1], np.log(trueval[2])], args=(x, y, yerr))
    m_ml, b_ml, lnf_ml = result["x"]



    ndim, nwalkers = 3, 100
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)] ## starting point

    pdb.set_trace()

    import emcee
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x,y,yerr))
    

    sampler.run_mcmc(pos,500)

    if showburn:
        ## to plot 'burn-in'
        ## samples = sampler.chain.reshape((-1,3))
        ax1 = subplot(311)
        ax2 = subplot(312)
        ax3 = subplot(313)
        for i in range(len(sampler.chain)):
            ax1.plot(range(len(sampler.chain[i])),sampler.chain[i,:,0],'k-',alpha=0.3)
            ax2.plot(range(len(sampler.chain[i])),sampler.chain[i,:,1],'k-',alpha=0.3)
            ax3.plot(range(len(sampler.chain[i])),sampler.chain[i,:,2],'k-',alpha=0.3)

        ax1.set_ylabel('$m$')
        ax2.set_ylabel('$b$')
        ax3.set_ylabel('$f$')
        ## ax1.set_xlim(0,500)
        ## ax2.set_xlim(0,500)
        ## ax3.set_xlim(0,500)
        u.adjust_spines(ax1,['left'])
        u.adjust_spines(ax2,['left'])
        u.adjust_spines(ax3,['left','bottom'])
        ax3.set_xlabel('step number')
        show()
        clf()


    if showcorner:
        samples = sampler.chain[:,50:,:].reshape((-1,ndim))
        ## to plot likelihood contours
        import corner ## using the Foreman-Mackey cornerplot util
        fig=corner.corner(samples, labels=['$m$','$b$','$\ln\,f$'],
                         truths=[trueval[0],trueval[1],log(trueval[2])])
        show()
        clf()
    

    samples = sampler.chain[:,50:,:].reshape((-1,ndim))
    samples[:, 2] = np.exp(samples[:, 2]) ## changes ln f back to f
    ## gives back the 68% percentile range of values on each dimension
    m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                 zip(*np.percentile(samples, [16, 50, 84],
                                                    axis=0))) 
    
    ax = subplot(111)
    ax.errorbar(x, y, yerr=yerr, fmt='o', color='red')
    xx = arange(-1, 11, 0.1)
    ax.plot(xx, linef(xx,*(m_ml,b_ml)),'b--')
    ax.plot(xx, linef(xx,*(m_mcmc[0],b_mcmc[0])),'b-')


    show()
    

    
