#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit
rcParams['font.size']=14.0
rcParams['figure.figsize']=9,6

from kapteyn import kmpfit
#matplotlib.rcParams.update({'font.size': 14})

myir = '#CC6677'
## myuv = '#117733'
myuv = '#DDCC77' #0.5
bar1 = myir
bar2 = '#88CCEE'

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

ax = subplot(111)
ax2 = ax.twinx()

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
h1 = ax.errorbar(redshifts, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
                 label = 'Madau & Dickinson (2014)',
                 fmt='o',mec=myuv,
                 mfc=myuv,ms = 5, c=myuv,
                 zorder=1, alpha=0.5)


compendium1 = zeros((len(data),5),)
compendium1[:,0]=redshifts
compendium1[:,1]=d_red
compendium1[:,2]=csfh
compendium1[:,3]=ep_csfh
compendium1[:,4]=em_csfh

                   

p0 = [0.015, 2.9, 2.7, 7.0]
popt, pcov = curve_fit(func, redshifts, csfh, p0 = p0, sigma = ep_csfh)
perr = sqrt(diag(pcov))
print(popt, perr)

fobj = kmpfit.simplefit(model,p0, redshifts, csfh, err= ep_csfh)
dfdp= dfdp_m(fobj.params,redshifts)

confprob = 0.68
yjunk, upperband, lowerband = fobj.confidence_band(redshifts, dfdp, confprob, model)


## data = loadtxt('md14_UV.txt')
## redshifts1 = 0.5*(data[:,1]+data[:,0])
## d_red = 0.5*(data[:,1]-data[:,0])
## csfh = 10**(data[:,2])
## csfha = csfh
## ep_csfh = 10**(data[:,2]+data[:,3])
## em_csfh = 10**(data[:,2]-data[:,4])


## ax.errorbar(redshifts1, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
##             label = 'FUV Tracers',
##             fmt='o',mec=myuv,
##             mfc=myuv,ms = 5, c=myuv,
##             zorder=2, alpha=0.5)
## data = loadtxt('md14_IR.txt')
## redshifts2 = 0.5*(data[:,1]+data[:,0])
## d_red = 0.5*(data[:,1]-data[:,0])
## csfh = 10**(data[:,2])
## csfhb = csfh
## ep_csfh = 10**(data[:,2]+data[:,3])
## em_csfh = 10**(data[:,2]-data[:,4])
## ax.errorbar(redshifts2, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
##             label='IR Tracers',ms=5,
##             fmt='o',color=myir,
##             zorder=10, alpha=0.5)


ax2.plot(redshifts, rz.MD_sfr_2014(redshifts, fink=False), 'b--', label=r'$\dot{\rho}_{\star}(z)$ Madau & Dickinson (2014)')
ax2.plot(redshifts, rz.MD_sfr_2014(redshifts), 'b:', label = r'$\dot{\rho}_{\star}(z)$  Finkelstein et al. (2014)')
## ax.plot(redshifts, func(redshifts, *fobj.params), 'b-', label = 'Fit to Madau & Dickinson',zorder=10)
## ax.fill_between(redshifts, upperband, lowerband, color = bar2, alpha = 0.4)


zz  = linspace(0., 8.0, 100)
dust = loadtxt('md14_AUV.txt')
p_dust, pcov = curve_fit(func, 0.5*(dust[:,0]+dust[:,1]), dust[:,2], p0 = p0)
perr = sqrt(diag(pcov))
print(p_dust, perr)


data = loadtxt('driver18.txt')
redshifts = 0.5*(data[:,2]+data[:,1])
d_red = 0.5*(data[:,2]-data[:,1])
csfh = 10**(data[:,3])

nscale = 1+10**(0.4*func(redshifts,*p_dust))
csfh=csfh*nscale*0.7**3

ep_csfh = 10**(log10(csfh)+sqrt(data[:,5]**2+data[:,6]**2+data[:,7]**2))
em_csfh = 10**(log10(csfh)-sqrt(data[:,5]**2+data[:,6]**2+data[:,7]**2))

compendium2 = zeros((len(data),5),)
compendium2[:,0]=redshifts
compendium2[:,1]=d_red
compendium2[:,2]=csfh
compendium2[:,3]=ep_csfh
compendium2[:,4]=em_csfh

p0 = [0.015, 2.9, 2.7, 7.0]
popt, pcov = curve_fit(func, redshifts, csfh, p0 = p0, sigma = ep_csfh)
perr = sqrt(diag(pcov))
print(popt, perr)

h2 = ax.errorbar(redshifts, csfh, xerr=d_red, yerr=[csfh-em_csfh,ep_csfh-csfh],
                 label=r'Driver et al. (2018)$^{a}$',
                 fmt='o',color='#999933',ms=7, #0.3
                 zorder=4, alpha=0.9)

## fobj = kmpfit.simplefit(model,p0, redshifts, csfh, err= ep_csfh)
## dfdp= dfdp_m(fobj.params,redshifts)

## confprob = 0.68
## yjunk, upperband, lowerband = fobj.confidence_band(redshifts, dfdp, confprob, model)

## ax.plot(redshifts, func(redshifts, *fobj.params), '-', color=myir, label='Fit to Driver')
## ax.fill_between(redshifts, upperband, lowerband, color = bar1, alpha=0.4)

redshifts=arange(0,10,0.1)
p0 = [0.015, 1.5, 5.0, 6.1]
perr = [0.001, 0.1, 0.2, 0.2]
ax2.plot(redshifts, func(redshifts, *p0)/0.7, 'g--', label = r'$\dot{\rho}_{\star}(z)$ CC SN rates, Strolger et al. (2015)', alpha=0.3)
zz = linspace(min(redshifts), max(redshifts), 100)
cov = diag(array(perr)**2.)
ps = np.random.multivariate_normal(p0, cov, 10000)
ysample = asarray([func(zz, *pi)/0.7 for pi in ps])
lower = percentile(ysample, 15.9, axis=0)
upper = percentile(ysample, 84.1, axis=0)
#ax2.fill_between(zz, upper, lower, color=myuv, alpha=0.2)


## from Wu et al. FIR background galaxies
p0 = [0.0157, 2.51, 3.64, 5.46]
perr = [0.0004, 0.04, 0.05, 0.1]
ax2.plot(redshifts, func(redshifts, *p0)*0.7, 'r--', label = r'$\dot{\rho}_{\star}(z)$ FIR Background, Wu et al. (2018)', alpha=0.3)


data = loadtxt('bouwens2015.txt')
redshifts = data[:,0]
csfh = 10**(data[:,1])*0.7
ep_csfh = 10**(log10(csfh)+data[:,2])
em_csfh = 10**(log10(csfh)-data[:,3])


compendium3 = zeros((len(data),5),)
compendium3[:,0]=redshifts
compendium3[:,1]+=1.0
compendium3[:,2]=csfh
compendium3[:,3]=ep_csfh
compendium3[:,4]=em_csfh

h3 = ax.errorbar(redshifts, csfh, yerr=[csfh-em_csfh,ep_csfh-csfh],
                 label='Bouwens et al. (2015)',
                 fmt='v',mfc='None', mec='0.3',ms=7,
                 zorder=4, alpha=0.9)

data = loadtxt('Khusanova2019.txt')
redshifts = data[:,0]
csfh = 10**(data[:,1])*0.7
ep_csfh = 10**(log10(csfh)+data[:,2])
em_csfh = 10**(log10(csfh)-data[:,3])


compendium4 = zeros((len(data),5),)
compendium4[:,0]=redshifts
compendium4[:,1]+=1.0
compendium4[:,2]=csfh
compendium4[:,3]=ep_csfh
compendium4[:,4]=em_csfh

h3 = ax.errorbar(redshifts, csfh, yerr=[csfh-em_csfh,ep_csfh-csfh],
                 label='Khusanova et al. (2019)',
                 fmt='^',mfc='None', mec='0.3',ms=7,
                 zorder=4, alpha=0.9)


compendium = concatenate((compendium1, compendium2,), axis=0)
compendium = array(sorted(compendium, key=lambda x: x[0]))

popt, pcov = curve_fit(func, compendium[:,0], compendium[:,2], p0 = p0, sigma = compendium[:,3])
perr = sqrt(diag(pcov))
print(popt, perr)

fobj = kmpfit.simplefit(model,p0, compendium[:,0], compendium[:,2], err= compendium[:,3])
dfdp= dfdp_m(fobj.params,compendium[:,0])
## print(fobj.params)
      

confprob = 0.68
yjunk, upperband, lowerband = fobj.confidence_band(compendium[:,0], dfdp, confprob, model)


ax2.plot(compendium[:,0], func(compendium[:,0], *fobj.params), '-', color='blue', label='Our fit to all data')
ax2.fill_between(compendium[:,0], upperband, lowerband, color='blue', alpha=0.2)

## handles,labels = ax.get_legend_handles_labels()
## print(labels)
## handles =[handles[6], handles[7],
##           handles[2], handles[3], handles[5],
##           handles[0], handles[1], handles[4]]
## labels =[labels[6], labels[7],
##           labels[2], labels[3], labels[5],
##           labels[0], labels[1], labels[4]]

ax.set_yscale('log')
ax.set_ylim(3e-3,1.1)
ax.set_xlim(-0.1,6.1)

ax2.set_yscale('log')
ax2.set_ylim(3e-3,1.1)
ax2.set_xlim(-0.1,6.1)
ax.set_xlabel(' Redshift')
ax.set_ylabel(r'$\dot{\rho}_{\star}$($M_{\odot}\, yr^{-1}\, Mpc^{-3}\, h_{70}$)')
lg1=ax.legend(loc=1, frameon=False)
lg2=ax2.legend(loc=3, frameon=False)


u.allblack2(ax,lg1)
#u.adjust_spines(ax,['left','bottom'])
u.allblack2(ax2,lg2)
#u.adjust_spines(ax2,['right','bottom'])

savefig('figure_csfh_today.png',transparent=True)


