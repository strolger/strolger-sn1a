#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit

from kapteyn import kmpfit

def func(z, *p):
    A,B,C,D=p
    sfr=A*((1+z)**C)/(1+((1+z)/B)**D)
    return(sfr)

def model(p, z):
    A,B,C,D=p
    sfr=A*((1+z)**C)/(1+((1+z)/B)**D)
    return(sfr)
    
def dfdp_m(p,z):
    A, B, C, D = p
    term1 =((1+z)**C)/(1+((1+z)/B)**D)
    term2 = A*D*((1+z)**C)*(((1+z)/B)**D)/(B*(1+((1+z)/B)**D)**2)
    term3 = (A*((1+z)**C)*log(1+z))/(1+((1+z)/B)**D)
    term4 = (A*((1+z)**C)*(((1+z)/B)**D)*log((1+z)/B))/((1+((1+z)/B)**D)**2.)
    return ([term1, term2,term3,term4])


data = loadtxt('md14.txt')

redshifts = 0.5*(data[:,1]+data[:,0])

dummy = zeros((len(data),len(data[0])+1),)
dummy[:,1:] = data
dummy[:,0] = redshifts
data = array(sorted(dummy, key=lambda x: x[0]))

redshifts = data[:,0]
d_red = 0.5*(data[:,2]-data[:,1])
csfh = 10**(data[:,3])
ep_csfh = 10**(data[:,3]+data[:,4])
em_csfh = 10**(data[:,3]-data[:,5])

p0 = [0.015, 2.9, 2.7, 7.0]
popt, pcov = curve_fit(func, redshifts, csfh, p0 = p0, sigma = ep_csfh)
perr = sqrt(diag(pcov))
print(popt, perr)

fobj = kmpfit.simplefit(model,p0, redshifts, csfh, err= ep_csfh)
dfdp= dfdp_m(fobj.params,redshifts)

confprob = 0.68
yjunk, upperband, lowerband = fobj.confidence_band(redshifts, dfdp, confprob, model)


ax = subplot(111)
data = loadtxt('md14_UV.txt')
redshifts1 = 0.5*(data[:,1]+data[:,0])
d_red = 0.5*(data[:,1]-data[:,0])
csfh = 10**(data[:,2])
ep_csfh = 10**(data[:,2]+data[:,3])
em_csfh = 10**(data[:,2]-data[:,4])
ax.errorbar(redshifts1, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
            label = 'MD14 FUV Tracers',
            fmt='o',color='purple', zorder=10)
data = loadtxt('md14_IR.txt')
redshifts1 = 0.5*(data[:,1]+data[:,0])
d_red = 0.5*(data[:,1]-data[:,0])
csfh = 10**(data[:,2])
ep_csfh = 10**(data[:,2]+data[:,3])
em_csfh = 10**(data[:,2]-data[:,4])
ax.errorbar(redshifts1, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
            label='MD14 IR Tracers',
            fmt='o',color='red', zorder=10)


ax.plot(redshifts, rz.MD_sfr_2014(redshifts, fink=False), 'b--', label='Madau & Dickinson (2014)')
## ax.plot(redshifts, rz.MD_sfr_2014(redshifts), 'b:', label = 'Finkelstein et al. 2016 ')
ax.plot(redshifts, func(redshifts, *fobj.params), 'b-', label = 'Our fit')
ax.fill_between(redshifts, upperband, lowerband, alpha = 0.4)



data = loadtxt('driver18.txt')
redshifts = 0.5*(data[:,2]+data[:,1])
d_red = 0.5*(data[:,2]-data[:,1])
csfh = 10**(data[:,3])
ep_csfh = 10**(data[:,3]+sqrt(data[:,5]**2+data[:,6]**2+data[:,7]**2))
em_csfh = 10**(data[:,3]-sqrt(data[:,5]**2+data[:,6]**2+data[:,7]**2))

p0 = [0.015, 2.9, 2.7, 7.0]
popt, pcov = curve_fit(func, redshifts, csfh, p0 = p0, sigma = ep_csfh)
perr = sqrt(diag(pcov))
print(popt, perr)

ax = subplot(111)
ax.errorbar(redshifts, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
            label='Driver et al. (2018) data',
            fmt='o',color='k', zorder=10)



fobj = kmpfit.simplefit(model,p0, redshifts, csfh, err= ep_csfh)
dfdp= dfdp_m(fobj.params,redshifts)

confprob = 0.68
yjunk, upperband, lowerband = fobj.confidence_band(redshifts, dfdp, confprob, model)

ax.plot(redshifts, func(redshifts, *fobj.params), 'r-', label='Our fit')
ax.fill_between(redshifts, upperband, lowerband, alpha=0.4)


p0 = [0.015, 1.5, 5.0, 6.1]
perr = [0.001, 0.1, 0.2, 0.2]
ax.plot(redshifts, func(redshifts, *p0), 'g-', label = 'CC rates')

zz = linspace(min(redshifts), max(redshifts), 100)
cov = diag(array(perr)**2.)
ps = np.random.multivariate_normal(p0, cov, 10000)
ysample = asarray([func(zz, *pi) for pi in ps])
lower = percentile(ysample, 15.9, axis=0)
upper = percentile(ysample, 84.1, axis=0)
ax.fill_between(zz, upper, lower, color='green', alpha=0.2)



ax.set_yscale('log')
ax.set_ylim(3e-3,1.1)
ax.set_xlim(-0.1,6.1)
ax.set_xlabel(' Redshift')
ax.set_ylabel(r'$\dot{\rho}_{\star}$($M_{\odot}\, yr^{-1}\, Mpc^{-3}$)')
ax.legend(loc=2, ncol=2, frameon=False)

savefig('figure_csfh_today.png')


