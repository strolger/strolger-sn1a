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

ndim, nwalkers, nsteps = 3, 60, 268

bounds = [
    [-15.,15.],
    [0.,15.],
    [-15.,15.]
    ]
    

def get_sfhs(verbose=True, delete=False):
    import tarfile
    
    sfh_file = 'SFH_file.tgz'
    if not os.path.isfile(sfh_file) or delete:
        sfhf = glob.glob('../ALLSFH_new_z/gsznnpas/*.dat')
        sfhs = {}
        for sfh in sfhf:
            sfhv = int(os.path.basename(sfh).split('.')[0])
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


def get_sfhs_old(verbose=True, delete=False):
    
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
    

def match_sne_hosts(gxycat='None',snecat='candels_south_sneia.txt', verbose=False):
    from strolger_util import convertdegsex as conv
    if verbose: print('Matching events to host galaxies...')
    f = open(snecat,'r')
    lines = f.readlines()
    f.close()
    output=[]
    if verbose: print('#SN\tHost\toffset(as)\tredshift')
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
    if sum(dtd) == 0.0:
        return (-np.inf)
    
    dt = sum(diff(sfh_data[:,0]))/(len(sfh_data[:,0])-1)
    rate_fn = zeros((len(sfh_data),2),)
    tmp = convolve(sfh_data[:,1], dtd, 'full')*dt*scale*frac_ia ## now convolved result in forward time
    rate_fn[:,1]=tmp[:len(dtd)]

    rate_fn[:,0]=sfh_data[:,0]
    rate_fn[:,1]=rate_fn[:,1]
    rate = simps(rate_fn[:,1],x=rate_fn[:,0])
            
    return(rate)


def lnprior(p):
    m, w, k = p
    if bounds[0][0] < m < bounds[0][1] and bounds[1][0] < w < bounds[1][1] and bounds[2][0] < k < bounds[2][1]:
        return 0.0
    else:
        return -np.inf



def lnlike(p):
    LL1 = 0.0
    LL2 = 0.0
    ## tcp = 2.0 ## in years
    for k in sfhs.keys():
        r_gxy = rate_per_galaxy(sfhs[k], p0=p) ## in number per year
        if ((r_gxy == 0. ) | (not isfinite(r_gxy))):
            #LL1a = -np.inf ## in the case that the DTD is unphysical, return neg inf
            #LL2a = -np.inf
            return(-np.inf)
        else:
            N_expected_Ia_gxy = r_gxy * tcp[k]
            LL1a=N_expected_Ia_gxy
            if k in ia_host_codex[:,0]:
                LL2a=log(N_expected_Ia_gxy)
            else:
                LL2a=0.0
                
    LL1+=LL1a; LL2+=LL2a
    LL = - LL1 + LL2
    if not isfinite(LL):
        return(-np.inf)
    return(LL)


def lnlike_old(p):
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
    p0 = (3.3, 0.5, 2.2)
    ## p0 = (-1.5, 0.5, -12.)
    timestamp=datetime.datetime.now().strftime('%Y%m%d%H%M')
    out_sampler = 'mc_sfh_%s.pkl' %timestamp
    ## out_sampler = 'mc_sfh_201806121453.pkl'
    verbose=True
    delete=False
    
    t0 = time.time()
    if verbose: print ('Loading SFHs...')
    sfhs = get_sfhs(verbose=verbose, delete=delete)

    candels_cat = loadtxt('../ALLSFH_new_z/CANDELS_GDSS_znew_avgal_radec.dat')

    ia_host_codex=match_sne_hosts(gxycat=candels_cat,verbose=verbose)

    if verbose: print ('Getting Redshifts...')
    redshifts={}
    tcp={}
    for item in candels_cat:
        redshifts[int(item[0])]=item[1]
        
    if verbose: print ('Getting control times...')
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
    in_list = list(redshifts.values())
    out_list = u.recast(in_list, 0., rds, yds)[1]
    for item in range(len(in_list)):
        tcp[item+1]=out_list[item]

    #lnlike(p0)
    #pdb.set_trace()
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
    
    
    
