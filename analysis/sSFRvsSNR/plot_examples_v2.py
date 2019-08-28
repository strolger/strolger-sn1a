#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
import  plot_sfr_snr as pss
from strolger_util import imf
import warnings
warnings.simplefilter('ignore',RuntimeWarning)
np.seterr(divide='ignore', invalid='ignore') 
from scipy.integrate import simps,quad
from scipy.optimize import curve_fit



def func1(x, *p):
    a,d,tau = p
    return((x+a)**d*exp(-(x+a)/tau))

def func2(x,*p):
    a,B,C,tau=p
    return((((a-x)/tau)**B + ((a-x)/tau)**(-C))**(-1))


def rate_per_galaxy(sfh_data, lbu=13.65, lbl=0.05, p0 = None,
                    frac_ia = 0.05,
                    ):
    
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    
    ii = where((sfh_data[:,0]>lbl)&(sfh_data[:,0]<=lbu)) # cuts the data range
    sfh_data = sfh_data[ii]

    sfh_data[:,0] = lbu-sfh_data[:,0][::-1]
    sfh_data[:,1] = sfh_data[:,1][::-1] ## now in forward time
    warnings.simplefilter('ignore',RuntimeWarning)
    dtd = rz.powerdtd(sfh_data[:,0], *p0)
    if sum(dtd) == 0.0:
        return (-np.inf)
    
    dt = sum(diff(sfh_data[:,0]))/(len(sfh_data[:,0])-1)
    rate_fn = zeros((len(sfh_data),2),)
    tmp = convolve(sfh_data[:,1], dtd, 'full')*dt*scale*frac_ia ## now convolved result in forward time
    rate_fn[:,1]=tmp[:len(dtd)]

    rate_fn[:,0]=sfh_data[:,0]
    rate_fn[:,1]=rate_fn[:,1]
    rate = rate_fn[-1,1] ## only need the current value
    
    ## return(rate_fn[-1,1])
    return(rate_fn)


if __name__=='__main__':

    ing = []
    lbu = 13.65
    xx = arange(0.05, lbu, 0.005)
    tt = lbu - xx
    p0a = (1., 4.7, 0.28)
    p0b = (lbu-4., 2.6, 0.8, 0.7)
    p0c = (-3, 4.7, 0.28)
    p0d = (lbu-8,2.6, 0.8, 0.7)
    ax = subplot(111)

    ax2 = axes([0.65, 0.2, 0.23, 0.2])
    ax2.set_title('SFH')
    ax2.set_xlabel('Lookback time (Gyr)')
    ax2.set_ylabel(r' M$_{\odot}$ yr$^{-1}$')

    data = zeros((len(xx), 2),)
    data[:,0] = xx ## lbt
    data[:,1] = 300*func1(xx, *p0a) ## an expression of SFH
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9*0.4 ## time cumulated mass in lookback order
    rdata = rate_per_galaxy(data, p0 = (-1.0, 1.0), frac_ia = 0.06) ## lft rates per galaxy
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1]) ## recasting in look back time

    ax3 = axes([0.23, 0.60, 0.23, 0.2])
    ax3.plot(data[:,0], mdata, 'C0-')
    ax3.set_title('cumulative mass')
    ax3.set_ylabel(r'M$_{\odot}$')
    ax3.set_xlabel('Lookback time (Gyr)')
    ax2.plot(data[:,0], data[:,1], 'C0-', label='"Grean Pea" Dwarf')
    #ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), color='C0',marker='_',s=1)#, c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)
    ax.plot(log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0], '*', ms=10, color='k', alpha=0.3)
    ing.append([log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0]])

    data = zeros((len(xx), 2),)
    data[:,0] = xx
    data[:,1] = 18*func2(xx, *p0b)
    data[:,1][isnan(data[:,1])]=0.0
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9
    rdata = rate_per_galaxy(data, p0 = (-1.0, 1.0), frac_ia = 0.06)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])

    ax3.plot(data[:,0], mdata, 'C1--')
    ax2.plot(data[:,0], data[:,1],'C1--', label='Ultra Diffuse Dwarf')
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), color='C1', marker='_', s=1)#, c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)
    ax.plot(log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0], '*', ms=10, color='k', alpha=0.3)
    ing.append([log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0]])

    data = zeros((len(xx), 2),)
    data[:,0] = xx
    ##data[:,1] = 2.25*func1(xx, *p0c)+0.18*func2(xx, *p0d)
    p_o = (7, 1.5, 1.5)
    data[:,1] = u.gauss(xx,*p_o) 
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9*0.8
    rdata = rate_per_galaxy(data, p0 = (-1.0, 1.0), frac_ia = 0.06)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])


    ax3.plot(data[:,0], mdata, 'C2:')
    ax2.plot(data[:,0],data[:,1],'C2:', label='Dwarf')
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), color='C2', marker='_', s=1)#,c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)
    ax.plot(log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0], '*', ms=10, color='k', alpha=0.3)
    ing.append([log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0]])


    data = zeros((len(xx), 2),)
    data[:,0] = xx
    ##data[:,1] = 2.25*func1(xx, *p0c)+0.18*func2(xx, *p0d)
    p_o = (7, 8., 4.5)
    data[:,1] = u.gauss(xx,*p_o) 
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9*0.8
    rdata = rate_per_galaxy(data, p0 = (-1.0, 1.0), frac_ia = 0.06)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])


    ax3.plot(data[:,0], mdata, 'C3:')
    ax2.plot(data[:,0],data[:,1],'C3:', label='Dwarf')
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), color='C3', marker='_',s=1)#c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)
    ax.plot(log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0], '*', ms=10, color='k', alpha=0.3)
    ing.append([log10(data[:,1]/mdata)[0], log10(rdn/mdata)[0]])



    ssfx = arange(-13, -7.5, 0.5)
    p_n = (1.19e-7, 0.586, 1.01e-11, 1.04e-9)
    ax.plot(ssfx, log10(pss.pieces(10**ssfx, *p_n)), 'k--', label='Andersen & Hjorth (2018) model')


    ## p_m = (0.12212169, 3.08357917, 5.71209495)
    ## p_e = (0.00970476, 0.19557979, 0.97953782)
    ## ##p_c = diag(array(p_e)**2.)
    ## p_c = array([[  9.41823879e-05 ,  1.89562431e-03,   9.46250054e-03],
    ##              [  1.89562431e-03,   3.82514562e-02,   1.91379154e-01],
    ##              [  9.46250054e-03,   1.91379154e-01,   9.59494348e-01]])


    p_m = (0.16241173, 3.82133413, 9.13755898)
    p_c = array([[4.38667393e-05, 8.82913042e-04, 4.40728962e-03],
                 [8.82913042e-04, 1.78161410e-02, 8.91374714e-02],
                 [4.40728962e-03, 8.91374714e-02, 4.46897689e-01]])
    
    ## ing = array(ing)
    ## p0 = (1.,1.,0.0)
    ## popt, pcov = curve_fit(pss.line_fn, ing[:,0], ing[:,1], p0=p0)


    ax.plot(ssfx, pss.line_fn(ssfx, *p_m), 'b-', label='Our model')
    ps = np.random.multivariate_normal(p_m, p_c, 10000)
    ysample = asarray([pss.line_fn(ssfx, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(ssfx, upper, lower, color='blue', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(ssfx, upper, lower, color='blue', alpha=0.2)
    

    ax.set_xlim(-12.5, -5.)
    ax.set_ylim(-14, -10.35)
    ax3.set_ylim(1,2e10)
    ax3.set_yscale('log')

    ax.set_ylabel('log(sSNR)')
    ax.set_xlabel('log(sSFR)')

    ax4 = axes([0.65, 0.6, 0.23, 0.2])
    p0 = (-1.0, 1.0)
    ax4.plot(xx, rz.powerdtd(xx, *p0))
    ax4.set_title('DTD')
    ax4.set_xlabel('Delay time (Gyr)')
    


    savefig('temp2.png')

    
