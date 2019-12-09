#!/usr/bin/env python
'''
This is the MCMC for SFHs and DTD
L. Strolger
5/2018
'''
import os,sys,pdb,scipy,glob,pickle,datetime
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
import control_time as tc
import warnings

## ndim, nwalkers, nsteps = 4, 30, 100
ndim, nwalkers, nsteps = 4, 100, 225

bounds = [
    #[0., 1.],
    [-10.,0.],
    [-2000., 2000.],
    [0.01,100.],
    [-500.,500.]
    ]


def get_sfhs(verbose=True, delete=False):
    import tarfile
    
    sfh_file = 'SFH_file.tgz'
    if not os.path.isfile(sfh_file) or delete:
        sfhf = glob.glob('../ALLSFH_new_z/g[sn]znnpas/*.dat')
        sfhs = {}
        for sfh in sfhf:
            if 'gsz' in sfh:
                sfhv = int(os.path.basename(sfh).split('.')[0])+200000
            else:
                sfhv = int(os.path.basename(sfh).split('.')[0])+100000
            sfhs[sfhv]=loadtxt(sfh)
        pickle.dump(sfhs, open(sfh_file.replace('tgz','pkl'),'wb'))
        tar = tarfile.open(sfh_file,mode='w:gz')
        tar.add(sfh_file.replace('tgz','pkl'))
        tar.close()
    else:
        tar=tarfile.open(sfh_file, mode='r:gz')
        tar.extractall()
        tar.close()
        sfhs=pickle.load(open(sfh_file.replace('tgz','pkl'),'rb'))
    os.remove(sfh_file.replace('tgz','pkl'))
    return(sfhs)


def match_sne_hosts(gxycat='None', verbose=False):
    from strolger_util import convertdegsex as conv
    if verbose: print('Matching events to host galaxies...')

    output=[]
    snecats=['candels_south_sneia.txt','candels_north_sneia.txt']
    if verbose: print('#SN\tHost\toffset(as)\tredshift')
    for snecat in snecats:
        f = open(snecat,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if line.startswith('#'): continue
            ra=line.split()[1]
            dec=line.split()[2]
            rr,dd = conv.s2d(ra,dec)
            offsets = sqrt((rr-gxycat[:,-2])**2+(dd-gxycat[:,-1])**2)
            idx = where(offsets == min(offsets))
            if min(offsets)*3600.0 < 10.0:
                if verbose: print('%s\t%d\t%2.2f\t%.1f' %(line.split()[0].strip('\t'),gxycat[idx][:,0],min(offsets)*3600.0,gxycat[idx][:,1]))
                output.append([gxycat[idx][:,0],line.split()[0].strip('\t'),min(offsets)*3600.0])
            else:
                if verbose: print('%s has no match' %(line.split()[0].strip('\t')))
    return(array(output))
    


def rate_per_galaxy(sfh_data, lbu=13.65, lbl=0.05, p0 = None,
                    ##frac_ia = 0.05,
                    ):
    
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale * 0.70**2.
    if not tuple(p0):
        p0 = (-3, -1.4, 3.5, -1.0)
    
    ii = where((sfh_data[:,0]>lbl)&(sfh_data[:,0]<=lbu)) # cuts the data range
    sfh_data = sfh_data[ii]

    sfh_data[:,0] = lbu-sfh_data[:,0][::-1]
    sfh_data[:,1] = sfh_data[:,1][::-1] ## now in forward time
    warnings.simplefilter('ignore',RuntimeWarning)
    dtd = rz.dtdfunc(sfh_data[:,0], *p0[1:])
    if sum(dtd) == 0.0:
        return (-np.inf)
    
    dt = sum(diff(sfh_data[:,0]))/(len(sfh_data[:,0])-1)
    rate_fn = zeros((len(sfh_data),2),)
    tmp = convolve(sfh_data[:,1], dtd, 'full')*dt*scale*exp(p0[0]) ## now convolved result in forward time
    rate_fn[:,1]=tmp[:len(dtd)]

    rate_fn[:,0]=sfh_data[:,0]
    rate_fn[:,1]=rate_fn[:,1]
    rate = rate_fn[-1,1] ## only need the current value
    
    return(rate_fn[-1,1])


def lnprior(p):
    f, m, w, k = p
    if bounds[0][0] < f < bounds[0][1] and bounds[1][0] < m < bounds[1][1] and bounds[2][0] < w < bounds[2][1] and bounds[3][0] < k < bounds[3][1]:
        return 0.0
    else:
        return -np.inf



def lnlike(p):
    LL1 = 0.0
    LL2 = 0.0
    for k in sfhs.keys():
        LL1a = 0.0
        LL2a = 0.0
        r_gxy = rate_per_galaxy(sfhs[k], p0=p)#, frac_ia=p[0]) ## in number per year
        if ((r_gxy == 0. ) | (not isfinite(r_gxy))):return(-np.inf)
        
        N_expected_Ia_gxy = r_gxy * tcp[k]
        LL1a=N_expected_Ia_gxy
        if k in ia_host_codex[:,0]:
            LL2a=log(N_expected_Ia_gxy)
        LL1+=LL1a; LL2+=LL2a
    LL = - LL1 + LL2
    if not isfinite(LL):
        return(-np.inf)
    return(LL)


def lnprob(p):
    lp = lnprior(p)
    if not isfinite(lp):
        return(-np.inf)
    else:
        return(lp + lnlike(p))


def choose_inbounds(p0, nwalkers, bounds):
    i = 0
    p_out=[]
    while i < nwalkers:
        ip0 = p0 + np.random.randn(len(p0))
        if isfinite(lnprior(ip0)):
            p_out.append(ip0)
            i+=1
    return(p_out)

if __name__ == '__main__':


    import time
    import multiprocessing as mpc

    ncore = mpc.cpu_count()
    step_sep = 1.0
    #p0 = (0.05, -1200., 50, 200)
    p0 = (-3, -1200., 50, 200)
    timestamp=datetime.datetime.now().strftime('%Y%m%d%H%M')
    out_sampler = 'mc_sfh_%s.pkl' %timestamp
    verbose=True
    delete=False
    
    t0 = time.time()
    if verbose: print ('Loading SFHs...')
    sfhs = get_sfhs(verbose=verbose, delete=delete)

    candels_cat_north = loadtxt('../ALLSFH_new_z/CANDELS_GDSN_znew_avgal_radec.dat')
    candels_cat_north = np.delete(candels_cat_north,[40,41],1) # removes two flag columns
    candels_cat_north[:,0]+=100000 #disambiguate the indexes
    candels_cat_south = loadtxt('../ALLSFH_new_z/CANDELS_GDSS_znew_avgal_radec.dat')
    candels_cat_south[:,0]+=200000 #disambiguate the indexes
    candels_cat = concatenate((candels_cat_north, candels_cat_south), axis=0)

    ia_host_codex=match_sne_hosts(gxycat=candels_cat,verbose=verbose)
    if verbose: print ('Getting Redshifts...')
    redshifts={}
    tcp={}
    for item in candels_cat:
        redshifts[int(item[0])]=item[1]
        
    if verbose: print ('Getting control times...')
    if not os.path.isfile('tcp.pkl') or delete:
        rds = arange(0.001, 5.5,0.1)
        yds = []
        for item in rds:
            tmp1 = tc.run(item,45.0,26.2,type=['ia'],dstep=3,dmstep=0.5,dastep=0.5,
                         verbose=False,plot=False,parallel=False,Nproc=1,
                         prev=0.0, extinction=False)*(1.0+item)
            tmp2 = tc.run(item,45.0,26.2,type=['ia'],dstep=3,dmstep=0.5,dastep=0.5,
                         verbose=False,plot=False,parallel=False,Nproc=1,
                         prev=45.0, extinction=False)*(1.0+item)
            tmp = 2*tmp1 + 8*tmp2
            yds.append(tmp)
        yds=array(yds)
        for k,v in redshifts.items():
            out_v = float(u.recast(v, 0., rds, yds)[1])
            tcp[k]=out_v
        pickle.dump(tcp,open('tcp.pkl','wb'))
    else:
        tcp = pickle.load(open('tcp.pkl','rb'))

    ## pdb.set_trace()
    if not os.path.isfile(out_sampler) or delete:
        if verbose: print ('Running MCMC... watch the load sparkle!')
        ## p0 = [p0 + step_sep*np.random.randn(ndim) for i in range(nwalkers)]
        p0 = choose_inbounds(p0,nwalkers,bounds)
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=ncore-1)
        warnings.simplefilter('ignore',RuntimeWarning)
        warnings.simplefilter('ignore',ValueError)
        sampler.run_mcmc(p0,nsteps)
        samples = sampler.chain
        pickle.dump(samples,open(out_sampler,'wb'))
    else:
        if verbose: print ('Reading MCMC results from %s...' %out_sampler)
        samples = pickle.load(open(out_sampler,'rb'))
    t1 = time.time()
    print(t1-t0)


    ## samples = samples[:,50:,:].reshape((-1,ndim))
    samples = samples.reshape((-1,ndim))
    ## gives back the 68% percentile range of values on each dimension
    f_mcmc, m_mcmc, w_mcmc, k_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                         zip(*np.percentile(samples, [16, 50, 84],
                                                            axis=0)))     
    print(r'parameters: $\varepsilon=%2.2f\pm%2.2f$; $\xi=%2.2f\pm%2.2f$; $\omega=%2.2f\pm%2.2f$; $\alpha=%2.2f\pm%2.2f$' %(f_mcmc[0], f_mcmc[1],
                                                                                                                m_mcmc[0], m_mcmc[1],
                                                                                                                w_mcmc[0], w_mcmc[1],
                                                                                                                k_mcmc[0], k_mcmc[1]))
    import corner
    samples = samples.reshape((-1,ndim))
    fig = corner.corner(samples,labels=[r'$\varepsilon$',r'$\xi$',r'$\omega$',r'$\alpha$'],
                        truths=[f_mcmc[0], m_mcmc[0], w_mcmc[0], k_mcmc[0]])
    fig.savefig('figure_preliminary_sfh_corners.png')
    
    
    
