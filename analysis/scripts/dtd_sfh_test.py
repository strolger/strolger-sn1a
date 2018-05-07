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
import warnings

def rate_per_galaxy(sfh, lbu=13.65, lbl=0.05, p0 = None,
        frac_ia = 0.05, plotit=False, title=None, testing=False,
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

    warnings.simplefilter('ignore',RuntimeWarning)
    dtd = rz.dtdfunc(sfh_data[:,0], *p0)

    if testing:
        ## makes the dtd function gaussian
        p1 = (1.,4.0, 0.01)
        dtd = u.gauss(sfh_data[:,0], *p1)
    

    rate_fn = zeros((len(sfh_data),2),)
    tmp = convolve(sfh_data[:,1], dtd, 'full')*scale*frac_ia ## now convolved result in forward time
    rate_fn[:,1]=tmp[:len(dtd)]

    rate_fn[:,0]=sfh_data[:,0]
    rate_fn[:,1]=rate_fn[:,1]
    rate = simps(rate_fn[:,1],x=rate_fn[:,0])

    if testing:
        ## a check of the shift in the effective time on each
        tmp1 = sum((sfh_data[:,0])*sfh_data[:,1])/sum(sfh_data[:,1])
        tmp2 = sum(rate_fn[:,1]*(rate_fn[:,0]))/sum(rate_fn[:,1])
        print(tmp1,tmp2,tmp2-tmp1) 

    if plotit:
        clf()
        ax1 = subplot(111)
        ax1.plot(lbu-sfh_data[:,0], sfh_data[:,1], '--', color='0.25')#,alpha=0.3)
        ax1.set_xlabel('Time (Gyr)')
        ax1.set_ylabel('$\psi(M_{\odot}\,yr^{-1})$')
        if title: ax1.set_title(title)

        ax3 = ax1.twinx()
        ax3.plot(lbu-rate_fn[:,0],rate_fn[:,1], 'k-', label = '$R=%2.2f\, yr^{-1}$' %(rate), alpha=0.3)
        ax3.set_ylabel('$\#\,yr^{-1}$')
        ax3.legend(frameon=False)

        ttt,ddd = rz.greggio()
        ttt = array(ttt); ddd = array(ddd)
        ii = where(ttt>0)
        ttt=ttt[ii]
        ddd=ddd[ii]
        time = arange(0,15,0.1)
        dtd = rz.dtdfunc(time,*p0)

        ax2 = axes([0.65, 0.6, 0.23, 0.2])
        ax2.plot(time,dtd,'b-', label= 'Fit')#label='Norm = %1.1f' %(simps(dtd,x=time)))
        ax2.plot(ttt,ddd,'b--', label='Greggio')
        ax2.set_ylabel('$\Phi$')
        ax2.set_xlabel('Delay Time (Gyr)')
        ax2.set_xlim(0,12)
        ax2.set_ylim(0,1.3)
        ax2.legend(frameon=False)

        if title:
            savefig(title+'.png')
        else:
            savefig('junk.png')
            
    return(rate)



def liklihood(N):
    LL1 = N
    LL2 = log(N)
    return(LL1,LL2)

if __name__=='__main__':


    p0 = (-2.3, 4.1, -11.2)
    p0 = (-4.3, 7.2, -4.0)
    ## p0 = None

    ## sfh  = '00529.dat'
    ## sfh  = '../SN_SFH/SFH_goodsn/'+sfh
    ## r_gxy = rate_per_galaxy(sfh, p0=p0, plotit=True) ## in number per year
    ## print(r_gxy)
    ## sys.exit()

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



    f = open('goods_1as.csv')
    lines = f.readlines()
    f.close()
    ias = []
    for line in lines:
        if line.startswith('#'): continue
        ias.append(line.split(',')[0])

        
    if os.path.isfile('output.txt'): os.remove('output.txt')
    f = open('output.txt','w')
    LL= 0.0
    for event in codex:
        index = '%05d' %int(event[0])
        try:
            sfh = glob.glob('../SN_SFH/*/%s.dat'%index)[0]
        except:
            print('%s not found, skipping...' %index)
            continue
        tcp = 10.0 ## in years
        if event[1] in ias:
            r_gxy = rate_per_galaxy(sfh,title=event[1], p0=p0) ## in number per year
            N_expected_Ia_gxy = r_gxy * tcp
            LL1=-N_expected_Ia_gxy
            LL2=log(N_expected_Ia_gxy)
        else:
            r_gxy = rate_per_galaxy(sfh,title=event[0], p0=p0) ## in number per year
            N_expected_Ia_gxy = r_gxy * tcp
            LL1=-N_expected_Ia_gxy
            LL2=0.0
        LL += LL1 + LL2
        f.write('%2.2f %2.2f %2.2f\n' %(LL1, LL2, LL))
    print ('Log Likelihood = %2.2f' %(LL))
    f.close()
    
    ## pdb.set_trace()
