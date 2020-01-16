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
import matplotlib.image as mpimg

import warnings

def rate_per_galaxy(sfh, lbu=13.65, lbl=0.05, p0 = None,
        frac_ia = 0.05, plotit=True, title=None, testing=False,
        ):
    
    scale = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale * 0.70**2.
    if not p0:
        p0 = (-1.4, 3.5, -1.0)
    
    sfh_data = loadtxt(sfh)
    ii = where((sfh_data[:,0]>lbl)&(sfh_data[:,0]<=lbu)) # cuts the data range
    sfh_data = sfh_data[ii]

    sfh_data[:,0] = lbu-sfh_data[:,0][::-1]
    sfh_data[:,1] = sfh_data[:,1][::-1] ## now in forward time
    sfh_data[:,2] = sfh_data[:,2][::-1] ## now in forward time

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
        ax1.plot(lbu-sfh_data[:,0], sfh_data[:,1], '--', color='0.25',alpha=0.3)
        #ax1.plot(lbu-sfh_data[:,0], sfh_data[:,1]+sfh_data[:,2], '--', color='r',alpha=0.3)
        #ax1.plot(lbu-sfh_data[:,0], sfh_data[:,1]-sfh_data[:,2], '--', color='r',alpha=0.3)
        ax1.fill_between(lbu-sfh_data[:,0], sfh_data[:,1]-sfh_data[:,2], sfh_data[:,1]+sfh_data[:,2], color='0.25', alpha=0.3)
        ax1.set_xlabel('Lookback Time (Gyr)')
        ax1.set_ylabel('$\psi(M_{\odot}\,yr^{-1})$')
        if title: ax1.set_title(title)
        ## ax1.axvline(x=0, color='r')

        ax3 = ax1.twinx()
        ax3.plot(lbu-rate_fn[:,0],rate_fn[:,1]*1.e3, 'g-', lw=3,  label = '$R_{Ia}(0)=%2.2f\, (1000\, yr)^{-1}$' %(rate*1.e3), alpha=0.3)
        ax3.set_ylabel('$R_{Ia}[\#\,(1000\, yr)^{-1}]$')
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
        ax2.plot(time,dtd,'r-', label= 'Best Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
        ax2.plot(time,rz.powerdtd(time, *pwrl), 'b--', label=r'$t^{%.1f}$'%(pwrl[0]))
        ax2.set_ylabel('$\Phi$')
        ax2.set_xlabel('Delay Time (Gyr)')
        ax2.set_xlim(0,12)
        ax2.set_ylim(0,1.3)
        ax2.legend(frameon=False)

        ax4 = axes([0.6, 0.2, 0.25, 0.25])
        img = mpimg.imread('107.jpg') ## Thoth
        ax4.imshow(img[75:135,50:140])
        #img = mpimg.imread('206.jpg') ## Borg
        #ax4.imshow(img[75:135,70:160])

        ax4.set_yticks([])
        ax4.set_xticks([])
        ax1.set_xlim(0,8)
        ax3.set_ylim(-1.3,20.5)
        
        if title:
            savefig(title+'.png')
        else:
            savefig('figure_sfh_demo_v1.pdf')
            
    return(rate)



def liklihood(N):
    LL1 = N
    LL2 = log(N)
    return(LL1,LL2)

if __name__=='__main__':


    p0 = (-1200., 50, 200)

    ## sfh0 = 106969 ## Borg
    sfh0 = 216531 ## Thoth
    if sfh0 < 200000:
        sfh = '../ALLSFH_new_z/gnznnpas/%05d.dat'%(sfh0 - 100000)
    else:
        sfh = '../ALLSFH_new_z/gsznnpas/%05d.dat'%(sfh0 - 200000)
        
    r_gxy = rate_per_galaxy(sfh, p0=p0, plotit=True, frac_ia=0.06) ## in number per year
    print('R_Ia = %2.2f per millennium' %(r_gxy*1e3))

    candels_cat_north = loadtxt('../ALLSFH_new_z/CANDELS_GDSN_znew_avgal_radec.dat')
    candels_cat_north = np.delete(candels_cat_north,[40,41],1) # removes two flag columns
    candels_cat_north[:,0]+=100000 #disambiguate the indexes
    candels_cat_south = loadtxt('../ALLSFH_new_z/CANDELS_GDSS_znew_avgal_radec.dat')
    candels_cat_south[:,0]+=200000 #disambiguate the indexes
    candels_cat = concatenate((candels_cat_north, candels_cat_south), axis=0)
    idx = where(candels_cat[:,0] == sfh0)

    redshift = candels_cat[:,1][idx]
    tmp1 = tc.run(redshift,45.0,26.2,type=['ia'],dstep=3,dmstep=0.5,dastep=0.5,
                  verbose=False,plot=False,parallel=False,Nproc=1,
                  prev=0.0, extinction=False)*(1.0+redshift)
    tmp2 = tc.run(redshift,45.0,26.2,type=['ia'],dstep=3,dmstep=0.5,dastep=0.5,
                  verbose=False,plot=False,parallel=False,Nproc=1,
                  prev=45.0, extinction=False)*(1.0+redshift)
    tcp = 2*tmp1 + 8*tmp2
    print('Control Time = %2.2f years in CANDELS' %tcp)

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
    print ('TCP= %f'%tcp)

    
    ## pdb.set_trace()
