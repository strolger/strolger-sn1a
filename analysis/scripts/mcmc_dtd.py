#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from astropy import convolution
from scipy import signal
from stsci.convolve import boxcar
from scipy.optimize import curve_fit
from scipy.integrate import simps,quad
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import imf
import emcee

import warnings#,exceptions
## warnings.simplefilter("error",RuntimeWarning)


def get_sfhs(verbose=True, delete=False):
    
    sfhfile = 'sfh_file.pkl'
    if not os.path.isfile(sfhfile) or delete:
        codex = []
        f = open('../SN_SFH/goodsn_match_idcand.dat')
        lines=f.readlines()
        f.close()
        for line in lines:
            if line.startswith('#'): continue
            c = line.split()
            codex.append(c[:2])
        f = open('../SN_SFH/goodss_match_idcand.dat')
        lines=f.readlines()
        f.close()
        for line in lines:
            if line.startswith('#'): continue
            c = line.split()
            codex.append(c[:2])
        codex=array(codex)

        sfhs={}
        f = open('goods_1as.csv')
        lines = f.readlines()
        f.close()

        ias = []
        for line in lines:
            if line.startswith('#'): continue
            ias.append(line.split(',')[0])
        for event in codex:
            index = '%05d' %int(event[0])
            try:
                sfh = glob.glob('../SN_SFH/*/%s.dat'%index)[0]
            except:
                if verbose: print('%s not found, skipping...' %index)
                continue

            sfhd = loadtxt(sfh)
            if event[1] in ias:
                try:
                    sfhs[1].append(sfhd)
                except:
                    sfhs[1]=[sfhd]
            else:
                try:
                    sfhs[0].append(sfhd)
                except:
                    sfhs[0]=[sfhd]
        pickle.dump(sfhs,open(sfhfile,'wb'))
    else:
        print('Loading %s' %(sfhfile))
        sfhs = pickle.load(open(sfhfile,'rb'))
    return(sfhs)
    



def rate_per_galaxy(sfh_data, lbu=13.65, lbl=0.05, p0 = None,
                    frac_ia = 0.05,
                    ):
    
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    if not tuple(p0):
        p0 = (-1.4, 3.5, -1.0)
    
    ii = where((sfh_data[:,0]>lbl)&(sfh_data[:,0]<=lbu)) # cuts the data range
    sfh_data = sfh_data[ii]

    sfh_data[:,0] = lbu-sfh_data[:,0][::-1]
    sfh_data[:,1] = sfh_data[:,1][::-1] ## now in forward time
    warnings.simplefilter('ignore',RuntimeWarning)
    dtd = rz.dtdfunc(sfh_data[:,0], *p0)
    rate_fn = zeros((len(sfh_data),2),)
    tmp = convolve(sfh_data[:,1], dtd, 'full')*scale*frac_ia ## now convolved result in forward time
    rate_fn[:,1]=tmp[:len(dtd)]

    rate_fn[:,0]=sfh_data[:,0]
    rate_fn[:,1]=rate_fn[:,1]
    rate = simps(rate_fn[:,1],x=rate_fn[:,0])
            
    return(rate)


def lnprior(p):
    m, w, k = p
    if -15.0 < m < 15 and 0.0 < w < 15.0 and -15.0 < k < 15.0:
        return 0.0
    else:
        return -np.inf



def lnlike(p):
    LL= 0.0
    tcp = 10.0 ## in years
    for k in sfhs.keys():
        for i in range(len(sfhs[k])):
            r_gxy = rate_per_galaxy(array(sfhs[k][i]), p0=p) ## in number per year
            N_expected_Ia_gxy = r_gxy * tcp
            LL1=-N_expected_Ia_gxy
            if k==1:
                LL2=log(N_expected_Ia_gxy)
            else:
                LL2=0.0
        LL+= LL1 + LL2
    if not isfinite(LL):
        return(-np.inf)
    return(LL)

def lnprob(p):
    lp = lnprior(p)
    if not isfinite(lp):
        return(-np.inf)
    else:
        return(lp + lnlike(p))


if __name__ == '__main__':


    import time
    import multiprocessing as mpc
    ncore = mpc.cpu_count()
    ndim, nwalkers, nsteps = 3, 100, 1500
    step_size = 1.0
    p0 = (3.5, 0.2, 2.2) #real starting point
    p0 = (-1.4, 3.5, -1.0)
    p0 = (-4.3, 7.2, -4.0)
    
    out_sampler = 'sampler.pkl'
    verbose=False
    delete=False
    
    ## backend not supported prior to emcee 3.0
    ## delete = True
    ## outfile = 'status.h5'
    ## if (os.path.isfile(outfile) and delete): os.remove(outfile)
    ## backend = emcee.backends.HDFBackend(outfile)
    t0 = time.time()
    sfhs = get_sfhs(verbose=verbose, delete=delete)
    ## backend.reset(nwalkers, ndim)
    if not os.path.isfile(out_sampler) or delete:
        p0 = [p0 + step_size*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=ncore-1)#,backend=backend)
        warnings.simplefilter('ignore',RuntimeWarning)
        warnings.simplefilter('ignore',ValueError)
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
    m_mcmc, w_mcmc, k_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                 zip(*np.percentile(samples, [16, 50, 84],
                                                    axis=0)))     
    print(r'parameters: $\xi=%2.2f\pm%2.2f$; $\omega=%2.2f\pm%2.2f$; $\alpha=%2.2f\pm%2.2f$' %(m_mcmc[0], m_mcmc[1]
                                                                                               ,w_mcmc[0], w_mcmc[1]
                                                                                               ,k_mcmc[0], k_mcmc[1]))
    import corner
    samples = samples.reshape((-1,ndim))
    fig = corner.corner(samples,labels=[r'$\xi$',r'$\omega$',r'$\alpha$'],
                        truths=[m_mcmc[0], w_mcmc[0], k_mcmc[0]])
    fig.savefig('temporary.png')
    
    
    
