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
    tmp = convolve(sfh, dtd, 'full')*ddt*p[0]#*frac*scale ## now convolved result in forward time
    rate_fn=tmp[:len(dtd)]
    ## print(p,sum(rate_fn))
    ## if sum(rate_fn)==0.:
    ##     pdb.set_trace()
    return(rate_fn)

    

if __name__=='__main__':

    dt = 0.05
    age=13.6
    tt = arange(dt,age,dt)
    p_val = [0.060, -1156., 58., 220.]
    p_err = [[0.003, -0.004],
             [839., -596.],
             [14., -16.],
             [191., -186.]
             ]

    ## fout, (ax, ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
    ax = subplot(111)
    ax2 = axes([0.55, 0.6, 0.33, 0.25])
    ax2.plot(tt, rz.dtdfunc(tt, *p_val[1:]), 'b-')

    p_err=sqrt(array(p_err)[:,1]**2.)
    ## p_err=sqrt(p_err[:,0]**2.+p_err[:,1]**2.)


    tta = linspace(min(tt), max(tt), 300)
    cov = diag(array(p_err)**2.)
    ## ps = np.random.multivariate_normal(p_val, cov, 10000)

    outfile = 'bfres0.pkl'
    if os.path.isfile(outfile):
        print ('%s exists' %outfile)
        ps=pickle.load(open(outfile,'rb'))
    else:
        ps = np.random.multivariate_normal(p_val, cov, 100)
        pickle.dump(ps, open(outfile,'wb'))
    print ('have ps')

    
    ## timestamp=datetime.datetime.now().strftime('%Y%m%d%H%M')
    ## outfile = 'bfres_%s.pkl' %timestamp
    outfile = 'bfres1.pkl'
    if os.path.isfile(outfile):
        print ('%s exists' %outfile)
        ysample=pickle.load(open(outfile,'rb'))
    else:
        ysample = asarray([rz.dtdfunc(tta, *pi[1:]) for pi in ps])
        pickle.dump(ysample, open(outfile,'wb'))
    print ('have ysample')

    ##idx = where((sum(ysample, axis=1)==0.)&(ysample[:,2]<=0.))
    idx = where(ysample[:,2]>0.)
    ysample = ysample[idx]

    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax2.fill_between(tta, upper, lower, color='green', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax2.fill_between(tta, upper, lower, color='green', alpha=0.2)

    pwrl=(-1.,1.)
    ax2.plot(tt,rz.powerdtd(tt, *pwrl), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))

    
    ax2.set_xlim(0,12)
    ax2.set_ylim(0,0.8)

    frac = 0.062
    
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

    outfile = 'bfres2.pkl'
    if os.path.isfile(outfile):
        print ('%s exists' %outfile)
        ysample=pickle.load(open(outfile,'rb'))
    else:
        ysample = asarray([run_model(tta, *pi) for pi in ps])*scale
        pickle.dump(ysample, open(outfile,'wb'))
    idx = where(ysample[:,2]>0.)
    ysample = ysample[idx]
        

    ax.plot(zz, run_model(tta, *p_val)*scale, 'b-')
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
    ax.plot(zz, jdud, 'b:')
    ## pdb.set_trace()

    ax.set_xlim(0,2.5)

    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'$10^{-4}$ SNe Ia per year per Mpc$^3$')
    ax2.set_ylabel('$\Phi$')
    ax2.set_xlabel('Delay Time (Gyr)')
    ax.set_title(r'$k=%.4f\,M_{\odot}^{-1},\,\,\varepsilon=%2.1f\%%$' %(scale_k, p_val[0]*100))
    savefig('figure_fit_demo_werr.png', transparent=True)



