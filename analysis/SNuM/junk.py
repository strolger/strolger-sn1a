#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.integrate import quad
from strolger_util import cosmotools as ct

rcParams['figure.figsize']=13,16
rcParams['font.size']=24.0

def csfh_mo(z, A, B, C, D):
    sfr=A*((1+z)**C)/(1+((1+z)/B)**D)
    hz=70.*(0.27*(1+z)**3+0.73)**(1/2.)
    
    return(sfr/(hz*(1+z)))
   
def dtdfit(time,*p):
    ff, aa, bb, cc = p
    scale = 0.0210
    scale = scale *0.7**2.*1e4
    par_model = [0.0134, 2.55, 3.3, 6.1]
    sfh = rz.csfh_time(time, *par_model)
    dt = sum(diff(time))/(len(time)-1)
    p1 = (aa, bb, cc)
    res = rz.dtdfunc(time, *p1,norm=True)
    tmp = convolve(sfh, res, 'full')
    return(ff*tmp[:len(time)]*dt*scale)

def dtdfit_pwrl(time,*p):
    ff, aa, bb= p
    scale = 0.0210
    scale = scale *0.7**2.*1e4
    par_model = [0.0134, 2.55, 3.3, 6.1]
    sfh = rz.csfh_time(time, *par_model)
    dt = sum(diff(time))/(len(time)-1)
    pwrl = (aa, bb)
    dud = rz.powerdtd(tt, *pwrl)
    tmp = convolve(sfh, dud, 'full')
    return(ff*tmp[:len(time)]*dt*scale)


if __name__=='__main__':


    tt = arange(0.05, 13.6, 0.05)
    lbt = 13.6 - tt
    zz = array([ct.cosmoz(x, ho=70) for x in lbt])
    #z = arange(0,10., 0.05)
    p0 = (0.0134, 2.55, 3.3, 6.1)
    ko = 0.0210
    
    R = 0.27
    h = 0.7
    
    out = []
    for ii in zz:
        out.append((1-R)*1e12*quad(csfh_mo, ii,inf, args=p0)[0])
    out=array(out)

    data = loadtxt('evol_cos_mass.txt')
    redz = 0.5*(data[:,1]+data[:,0])
    ax3 = subplot(212)
    ax2 = subplot(211)
    ## ax1 = subplot(311)
    ## ax1.errorbar(redz, pow(10, data[:,2]), xerr=(data[:,1]-redz), yerr=[pow(10,data[:,2]+data[:,4])-pow(10,data[:,2]),
    ##                                                                     pow(10,data[:,2])-pow(10,data[:,2]-data[:,4])],
    ##              fmt='o')
    ## ax1.plot(zz, out*h, 'r-', label='Evolution of Stellar Mass Density')


    data1 = loadtxt('snia.txt')
    ax2.errorbar(data1[:,0], data1[:,2], xerr=data1[:,1],
                 yerr=[data1[:,3],data1[:,4]],
                 fmt='o', mfc='C3', mec='C3', color='C3')
    pc = [0.062, -1518, 51, 50]
    snriaa = dtdfit(tt,*pc)
    pc = [0.062, -1, 1]
    snriab = dtdfit_pwrl(tt,*pc)
    ax2.plot(zz, snriab, 'r-', label=R'SN Ia Rate History, $\beta=-1$')
    ax2.plot(zz, snriaa, 'r--', label=R'SN Ia Rate History, exponential')

    data2 = loadtxt('ccsne.txt')
    p1 = [0.015, 1.5, 5.0, 6.1]
    k1 = 0.0091
    ax2.plot(zz, rz.csfh(zz,*p0)/k1*0.7**2., 'C0-', label='CCSN Rate History')
    ax2.plot(zz, rz.csfh(zz,*p1)/k1*0.7**2., 'C0--', label=r'Fit to R$_{CC}$ data')
    ax2.errorbar(data2[:,0], data2[:,2], xerr=data2[:,1],
                 yerr=[data2[:,4],data2[:,3]],
                 fmt='o')

    ax3.plot(zz, (snriaa*1e-4)/(out*h)*100*1e10, 'C3--', lw=3, alpha=0.6, zorder=9, label='SNe Ia, exponential')
    ax3.plot(zz, (snriab*1e-4)/(out*h)*100*1e10, 'C3-', alpha=1, zorder=8, lw=3, label=r'SNe Ia, $\beta=-1$')
    
    ax3.plot(zz, (rz.csfh(zz,*p1)/k1*0.7**2*1e-4)/(out*h)*100*1e10, 'C0--', alpha=0.6, lw=3, zorder=7, label='CCSN fit')
    ax3.plot(zz, (rz.csfh(zz,*p0)/k1*0.7**2*1e-4)/(out*h)*100*1e10, 'C0-', alpha=0.6, zorder=3, lw=3, label='CCSN expected')

    ## ax3.plot(zz,((snriaa*1e-4)+(rz.csfh(zz,*p1)/k1*0.7**2*1e-4))/(out*h)*100*1e10,'--',
    ##          color='purple', lw=3,zorder=5)

    ## ## ax3.fill_between(zz,((snriab*1e-4)+(rz.csfh(zz,*p1)/k1*0.7**2*1e-4))/(out*h)*100*1e10, color='C4', alpha=1., zorder=4)
    ax3.plot(zz,((snriaa*1e-4)+(rz.csfh(zz,*p1)/k1*0.7**2*1e-4))/(out*h)*100*1e10,'C2--',
             alpha=0.3, zorder=4, lw=3, label='Low all SN est.')
    ax3.plot(zz,((snriaa*1e-4)+(rz.csfh(zz,*p0)/k1*0.7**2*1e-4))/(out*h)*100*1e10,'C2-', lw=3, alpha=0.3, zorder=2,
             label=r'Expected SNe per century per 10$^{10}$M$_{\odot}$')




    ## ## ax3.fill_between(zz, (rz.csfh(zz,*p0)/k1*0.7**2*1e-4)/(out*h)*100*1e10, color='b', alpha=0.6, zorder=3)
    ## ## ax3.fill_between(zz, (rz.csfh(zz,*p1)/k1*0.7**2*1e-4)/(out*h)*100*1e10, color='C4', alpha=1, zorder=4)


    ## #x0 = [0.0]
    ## #junk, out = u.recast(x0, 0.0, zz, ((snriab*1e-4)+(rz.csfh(zz,*p0)/k1*0.7**2*1e-4))/(out*h)*100*1e10)
    ## #pdb.set_trace()
    
    ## ax3.plot(zz,((snriaa*1e-4)+(rz.csfh(zz,*p0)/k1*0.7**2*1e-4))/(out*h)*100*1e10, 'C2--', lw=3, alpha=1, zorder=2)
    ## ## ax3.fill_between(zz,((snriab*1e-4)+(rz.csfh(zz,*p0)/k1*0.7**2*1e-4))/(out*h)*100*1e10, color='C2', alpha=0.3, zorder=1)

    ax3.set_xlim(0,3)
    ax2.set_xlim(0,3)
    ## ax1.set_xlim(0,3)
    ax3.set_ylim(0.01,9)
    ax2.set_ylim(0.09,11)
    ## ax1.set_ylim(1.1e7,1e9)
    ax2.set_xticklabels([])
    ## ax1.set_xticklabels([])
    
    subplots_adjust(hspace=0)
    ## ax1.set_yscale('log')
    ax2.set_yscale('log')
    ## ax3.set_yscale('log')

    ## ax1.legend(loc=8,frameon=False)
    ax2.legend(loc=8,ncol=2,frameon=False,fontsize=16)
    handles, labels = ax3.get_legend_handles_labels()
    ax3.legend(handles[::-1], labels[::-1], loc=2, frameon=False, fontsize=20)
    ## ax3.legend(loc=2,frameon=False)
    ax3.set_xlabel('Redshift')
    ## ax1.set_ylabel(r'$\rho_{\star}$ (M$_{\odot}$ Mpc$^{-3}$ $h_{70}$)')
    ax2.set_ylabel(r'R (10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ $h_{70}^3$)')
    ax3.set_ylabel('SNuM')
    ## ax1.set_title ('Type Ia Supernovae')
    savefig('SNuM_all.pdf')
    
