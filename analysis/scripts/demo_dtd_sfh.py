#!/usr/bin/env python
'''
This code runs the DTD test for a single choice of dtd parameters
'''
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
import control_time as tc

import warnings

def rate_per_galaxy(sfh, lbu=13.65, lbl=0.05, p0 = None,
        frac_ia = 0.05, plotit=True, title=None, testing=False,
        ):
    
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    if not p0:
        p0 = (-1.4, 3.5, -1.0)
    
    sfh_data = loadtxt(sfh)
    ii = where((sfh_data[:,0]>lbl)&(sfh_data[:,0]<=lbu)) # cuts the data range
    sfh_data = sfh_data[ii]

    sfh_data[:,0] = lbu-sfh_data[:,0][::-1]
    sfh_data[:,1] = sfh_data[:,1][::-1] ## now in forward time

    if testing:
        ## reframes everything in a more continuous time.
        n_out = IUS(sfh_data[:,0], sfh_data[:,1])
        time = arange(0,lbu+0.1,0.1)
        sfh_data=zeros((len(time),2),)
        sfh_data[:,0]=time
        sfh_data[:,1]=n_out(time) ## now continuous, evenly spaced, not sure that's important.

    if testing:
        ## just makes a SFH box
        box_ll=4.
        ii = where((sfh_data[:,0]>box_ll)&(sfh_data[:,0]<box_ll+2.))
        sfh_data[:,1] = 0.0
        sfh_data[:,1][ii]=10.

    #warnings.simplefilter('ignore',RuntimeWarning)
    dtd = rz.dtdfunc(sfh_data[:,0], *p0)

    if testing:
        ## makes the dtd function gaussian
        p1 = (1.,4.0, 0.01)
        dtd = u.gauss(sfh_data[:,0], *p1)
    

    dt = sum(diff(sfh_data[:,0]))/(len(sfh_data[:,0])-1)
    rate_fn = zeros((len(sfh_data),2),)
    tmp = convolve(sfh_data[:,1], dtd, 'full')*dt*scale*frac_ia ## now convolved result in forward time
    rate_fn[:,1]=tmp[:len(dtd)]#*concatenate((array([0.]),diff(sfh_data[:,0])),)

    rate_fn[:,0]=sfh_data[:,0]
    rate_fn[:,1]=rate_fn[:,1]
    rate = rate_fn[-1,1]

    if testing:
        ## a check of the shift in the effective time on each
        tmp1 = sum((sfh_data[:,0])*sfh_data[:,1])/sum(sfh_data[:,1])
        tmp2 = sum(rate_fn[:,1]*(rate_fn[:,0]))/sum(rate_fn[:,1])
        print(tmp1,tmp2,tmp2-tmp1) 

    if plotit:
        clf()
        ax1 = subplot(111)
        ax1.plot(lbu-sfh_data[:,0], sfh_data[:,1], '--', color='0.25')#,alpha=0.3)
        ax1.set_xlabel('Lookback Time (Gyr)')
        ax1.set_ylabel('$\psi(M_{\odot}\,yr^{-1})$')
        if title: ax1.set_title(title)

        ax3 = ax1.twinx()
        ax3.plot(lbu-rate_fn[:,0],rate_fn[:,1]*1.e3, 'k-', label = '$R_{Ia}(0)=%2.2f\, (1000\, yr)^{-1}$' %(rate*1.e3), alpha=0.3)
        ax3.set_ylabel('$R_{Ia}(\#\,(1000\, yr)^{-1})$')
        ax3.legend(frameon=False)

        ttt,ddd = rz.greggio()
        ttt = array(ttt); ddd = array(ddd)
        ii = where(ttt>0)
        ttt=ttt[ii]
        ddd=ddd[ii]
        time = arange(0.1,15,0.1)
        dtd = rz.dtdfunc(time,*p0)
    
        pwrl = (-1.0,1.0)
        ax2 = axes([0.6, 0.55, 0.25, 0.25])
        ax2.plot(time,dtd,'r-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
        ax2.plot(time,rz.powerdtd(time, *pwrl), 'b:', label=r'$t^{%.1f}$'%(pwrl[0]))
        ax2.plot(ttt,ddd,'b--', label='Greggio')
        ax2.set_ylabel('$\Phi$')
        ax2.set_xlabel('Delay Time (Gyr)')
        ax2.set_xlim(0,12)
        ax2.set_ylim(0,1.3)
        ax2.legend(frameon=False)
        #tight_layout()

        if title:
            savefig(title+'.png')
        else:
            savefig('figure_sfh_demo.png')
            
    return(rate)



def liklihood(N):
    LL1 = N
    LL2 = log(N)
    return(LL1,LL2)

if __name__=='__main__':


    p0 = (1.03, 7.98, 1.2)
    p0 = (3.3, 0.5, 2.2)
    p0 = (-6.2, 9.7, -0.7)

    sfh0 = 19962
    sfh = '../ALLSFH_new_z/gsznnpas/'+str(sfh0)+'.dat'
    r_gxy = rate_per_galaxy(sfh, p0=p0, plotit=True) ## in number per year
    print('R_Ia = %2.2f per mileneum' %(r_gxy*1e3))

    candels_cat = loadtxt('../ALLSFH_new_z/CANDELS_GDSS_znew_avgal_radec.dat')
    redshift = candels_cat[sfh0-1,1]
    tmp1 = tc.run(redshift,45.0,26.2,type=['ia'],dstep=3,dmstep=0.5,dastep=0.5,
                  verbose=False,plot=False,parallel=False,Nproc=1,
                  prev=0.0, extinction=False)*(1.0+redshift)
    tmp2 = tc.run(redshift,45.0,26.2,type=['ia'],dstep=3,dmstep=0.5,dastep=0.5,
                  verbose=False,plot=False,parallel=False,Nproc=1,
                  prev=45.0, extinction=False)*(1.0+redshift)
    tcp = 2*tmp1 + 8*tmp2
    print('Control Time = %2.2f years' %tcp)

    iahost=True
    if iahost:
        N_expected_Ia_gxy = r_gxy * tcp
        LL1=-N_expected_Ia_gxy
        LL2=log(N_expected_Ia_gxy)
    else:
        N_expected_Ia_gxy = r_gxy * tcp
        LL1=-N_expected_Ia_gxy
        LL2=0.0
    LL = LL1 + LL2
    print ('lnL for galaxies is %2.1e, and for Hosts is %2.1f' %(LL1, LL2))
    print ('Log Likelihood = %2.2f' %(LL))
    
    ## pdb.set_trace()
