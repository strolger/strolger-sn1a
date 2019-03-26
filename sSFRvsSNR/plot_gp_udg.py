#!/usr/bin/env python
import os,sys,pdb,scipy,glob
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
import warnings
warnings.simplefilter('ignore',RuntimeWarning)
np.seterr(divide='ignore', invalid='ignore') 


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
    #scale = scale * 0.70**2.
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
    rate = rate_fn[-1,1] ## only need the current value
    
    ## return(rate_fn[-1,1])
    return(rate_fn)




if __name__=='__main__':



    
    lbu = 13.65
    xx = arange(0.05, lbu, 0.005)
    tt = lbu - xx
    p0a = (1., 4.7, 0.28)
    p0b = (lbu, 2.6, 0.8, 0.7)
    p0c = (-5, 4.7, 0.28)
    p0d = (lbu-5,2.6, 0.8, 0.7)
    ax = subplot(111)

    ax2 = axes([0.65, 0.2, 0.23, 0.2])
    ax2.plot(xx,300*func1(xx, *p0a),'r-', label='"Grean Pea" Dwarf')
    ax2.plot(xx,300*func1(xx, *p0c)+18*func2(xx, *p0d),'r:', label='Dwarf')
    ax2.plot(xx,18*func2(xx, *p0b),'r--', label='Ultra Diffuse Dwarf')
    ax2.set_title('SFH')
    ax2.set_xlabel('Lookback time (Gyr)')
    ax2.set_ylabel(r' M$_{\odot}$ yr$^{-1}$')
    #ax2.legend(frameon=False)

    data = zeros((len(xx), 2),)
    data[:,0] = xx
    data[:,1] = 300*func1(xx, *p0a)
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9
    mdata[where(mdata<1e-5)]=0.0
    rdata = rate_per_galaxy(data, p0 = (-1258, 59, 248), frac_ia = 0.06)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])

    ax3 = axes([0.25, 0.2, 0.23, 0.2])
    ax3.plot(data[:,0], log10(data[:,1]/mdata))
    ax3.set_title('cumulative mass')
    ax3.set_xlabel('Lookback time (Gyr)')
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)


    data = zeros((len(xx), 2),)
    data[:,0] = xx
    data[:,1] = 18*func2(xx, *p0b)
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9
    mdata[where(mdata<1e-5)]=0.0
    rdata = rate_per_galaxy(data, p0 = (-1258, 59, 248), frac_ia = 0.06)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])

    ax3.plot(data[:,0], log10(data[:,1]/mdata), '--')
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)

    data = zeros((len(xx), 2),)
    data[:,0] = xx
    data[:,1] = 300*func1(xx, *p0c)+18*func2(xx, *p0d)
    data[:,1][isnan(data[:,1])]=0.0
    mdata = cumsum(data[:,1][::-1])[::-1]*0.005*1e9
    mdata[where(mdata<1e-5)]=0.0
    rdata = rate_per_galaxy(data, p0 = (-1258, 59, 248), frac_ia = 0.06)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])

    ax3.plot(data[:,0], log10(data[:,1]/mdata), ':')
    ax.scatter(log10(data[:,1]/mdata), log10(rdn/mdata), c=range(len(data[:,0])), cmap='plasma', marker='_', s=1)



    #ax.set_xlim(-15.5, -8.)
    #ax.set_ylim(-15.6, -13.0)
    savefig('temp.png')

    
