#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle,datetime
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import imf
from strolger_util import cosmotools as ct
from copy import copy, deepcopy
from scipy.integrate import simps,quad
import emcee

'''
This is the MCMC for SFDH and DTD
L. Strolger
5/2018
'''

def dtdfit(time,*p):
    ff, m, w, k = p
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale *0.7**2.*1e4

    par_model = [0.0134, 2.55, 3.3, 6.1]  ## Compendium
    ## par_model = [0.013, 2.6, 3.2, 6.1] ## MD14
    ## par_model = [0.009, 2.7, 2.5, 4.1] ## Driver18
    
    sfh = rz.csfh_time(time, *par_model)
    dt = sum(diff(time))/(len(time)-1)
    p1 = (m, w, k)
    res = rz.dtdfunc(time, *p1, norm=True)
    tmp = convolve(sfh, res, 'full')
    return(np.exp(ff)*tmp[:len(time)]*dt*scale)


def lnprior(p):
    ff, m, w, k, lnf = p
    if -10.0 < ff < 0.0 and -2000.0 < m < 2000.0 and 0.01 < w < 100.0 and -500.0 < k < 500.0 and -4.0 < lnf < 0.0:
        return 0.0
    else:
        return -np.inf

def lnlike(p, x, y, yerr):
    ff, m, w, k, lnf = p
    p1 = (ff, m, w, k)
    model = dtdfit(x,*p1)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    if sum(model)==0.:
        retval=-np.inf
    else:
        retval =  -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
    return retval

def lnprob(p, x, y, yerr):
    lp = lnprior(p)
    if not isfinite(lp):
        return(-np.inf)
    else:
        return(lp + lnlike(p, x, y, yerr))




if __name__=='__main__':

    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,verbose=False,splits=arange(0,1.167,0.167).tolist())
    data = deepcopy(brates)
    tt = 13.6-array([ct.cosmotime(x) for x in data[:,0]])
    data[:,0] = tt
    data = data[argsort(data[:,0])]
    
    import time
    import multiprocessing as mpc
    ncore = mpc.cpu_count()
    ndim, nwalkers, nsteps = 5, 100, 2000
    step_size = 1.0
    ## p0 = (0.05, 3.5, 2.5, 2.5, 0.01)
    ## p0 = (0.06, -100., 50, 20, -2.5)
    p0 = (-2.81, -1200, 50, 200, -2.5)

    mckaps = glob.glob('mc_sfd_*.pkl')
    if mckaps:
        timestamps = [int(x.split('.')[0].split('_')[-1]) for x in mckaps]
        idx = where(array(timestamps)==max(array(timestamps)))
        out_sampler = array(mckaps)[idx][0]
    else:
        timestamp=datetime.datetime.now().strftime('%Y%m%d%H%M')
        out_sampler = 'mc_sfd_%s.pkl' %timestamp

    verbose=False
    delete=False

    t0 = time.time()
    if not os.path.isfile(out_sampler) or delete:
        p0 = [p0 + step_size*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=ncore-1, args=(data[:,0], data[:,1], data[:,3]))
        sampler.run_mcmc(p0,nsteps)
        samples = sampler.chain
        pickle.dump(samples,open(out_sampler,'wb'))
    else:
        samples = pickle.load(open(out_sampler,'rb'))
    t1 = time.time()
    print(t1-t0)

    samples = samples[:,50:,:].reshape((-1,ndim))
    ## samples = samples.reshape((-1,ndim))
    ## gives back the 68% percentile range of values on each dimension
    ff_mcmc, m_mcmc, w_mcmc, k_mcmc, lnf_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                                    zip(*np.percentile(samples, [16, 50, 84],
                                                                       axis=0)))     
    print(r'scale: $f=%2.2f\pm%2.2f$'%(ff_mcmc[0], ff_mcmc[1]))
    print(r'parameters: $\xi=%2.2f\pm%2.2f$; $\omega=%2.2f\pm%2.2f$; $\alpha=%2.2f\pm%2.2f$' %(m_mcmc[0], m_mcmc[1]
                                                                                               ,w_mcmc[0], w_mcmc[1]
                                                                                               ,k_mcmc[0], k_mcmc[1]))

    md0 = u.binmode(samples[:,1],bins=20)[0]
    md1 = u.binmode(samples[:,2],bins=20)[0]
    md2 = u.binmode(samples[:,3],bins=20)[0]
    print(md0, md1, md2)

    import corner
    fig = corner.corner(samples,labels=['frac', r'$\xi$',r'$\omega$',r'$\alpha$', 'lnf'],
                        truths=[ff_mcmc[0], md0, md1, md2, lnf_mcmc[0]])
                        ## truths=[ff_mcmc[0], m_mcmc[0], w_mcmc[0], k_mcmc[0], lnf_mcmc[0]])
    fig.savefig('figure_preliminary_sfd_corners.png')
