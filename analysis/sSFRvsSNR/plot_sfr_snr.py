#!/usr/bin/env python
import os,sys,pdb,scipy,glob, pickle
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
#rcParams['font.size']=15.0

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


def peices(sSFR, *p):
    a, b,  x1, x2 = p
    yy = a*sSFR**b
    yy[where(sSFR < x1)]=a*x1**b
    yy[where(sSFR >= x2)]=a*x2**b
    return(yy)

def line_fn(x, *p):
    a,b,c = p
    y = a*x**2 +b*x+c 
    return(y)



def get_data(data, p0=(-1.4, 1., 1.), lbu=13.65, frac_ia=0.05, plotit=False):
    rdata = rate_per_galaxy(data, p0=p0, frac_ia=frac_ia)
    junk, rdn=u.recast(data[:,0], 0., lbu-rdata[:,0], rdata[:,1])
    out = zeros((len(data), 3),)
    out[:,0] = data[:,0]
    out[:,1] = data[:,1]/data[:,3]
    out[:,2] = rdn/data[:,3]
    if plotit:
        clf()
        ax1 = subplot(411)
        ax1.plot(data[:,0], data[:,1], 'r-', label='SFR')
        ax2 = subplot(412)
        ax2.plot(data[:,0], data[:,3]/1e10, 'b-', label='Total Mass')
        ax3 = subplot(413)
        ssfr = data[:,1]/data[:,3]*1e10
        ax3.plot(data[:,0], data[:,1]/data[:,3]*1e10, 'g-', label='sSFR')
        ax4 = subplot(414)
        ax4.plot(lbu-rdata[:,0], rdata[:,1], '-', color='purple', label='sSNR')
        ax1.set_xlim(0,20)
        ax2.set_xlim(0,20)
        ax3.set_xlim(0,20)
        ax4.set_xlim(0,20)
        ax4.set_xlabel('Lookback Time(Gyr)')
        ax1.legend()
        ax2.legend()
        ax3.legend()
        ax4.legend()
        savefig('figure_%06d_history.png'%(idx))

    return(out)



if __name__=='__main__':

    if not os.path.isfile('ssSFRN.pkl'):
        idxs = loadtxt('../analysis/scripts/host_idxs.txt')
        sfhs = {}
        for idx in idxs:
            print('%d' %idx)
            if idx > 200000:
                file = '../analysis/ALLSFH_new_z/gsznnpas/%05d.dat'%(idx-200000)
            else:
                file = '../analysis/ALLSFH_new_z/gnznnpas/%05d.dat'%(idx-100000)
            if os.path.isfile(file):
                data = loadtxt(file)
            else:
                print('%s not found'%file)
                continue
            ndata = get_data(data, p0=(-1258, 59, 248), frac_ia = 0.06, plotit=True) ## my model
            ## ndata = get_data(data, p0=(-23.1, 4.1, -1.1), frac_ia = 0.06, plotit=True) ## t**-1
            ## ndata = get_data(data, p0=(-13.4, 2.4, -1.1), frac_ia = 0.06, plotit=True) ## t**-1.1
            ## ndata = get_data(data, p0=(-4.3, 1.3, -1.8), frac_ia = 0.06, plotit=True) ## t**-1.1
            ## ndata = get_data(data, p0=(2.5, 1.0, 2.), frac_ia = 0.06, plotit=True) ## very delayed
            
            sfhs[idx]=ndata
            pickle.dump(sfhs,open('ssSFRN.pkl','wb'))
    else:
        sfhs = pickle.load(open('ssSFRN.pkl','rb'))
        

    data1 = loadtxt('../analysis/ALLSFH_new_z/Cami_GOODS-S_zbest.dat')
    data2 = loadtxt('../analysis/ALLSFH_new_z/Cami_GOODS-N_zbest.dat')
    

    ax = subplot(111)
    k1 = list(sfhs.keys())[0]
    ii = max(where(sfhs[k1][:,0]<0.05)[0])
    out = []
    for iii, idx in enumerate(sfhs.keys()):
        xx = log10(sfhs[idx][ii:,1])
        yy = log10(sfhs[idx][ii:,2])
        y0 = sfhs[idx][0,2]
        ## if iii == 0:
        ##     j1 = ax.scatter(xx, yy, c=range(len(xx)), cmap='plasma', marker='_', s=1, label='Temporal track')
        ## else:
        ##     j1 = ax.scatter(xx, yy, c=range(len(xx)), cmap='plasma', marker='_', s=1)
            

        if idx > 200000:
            x0 = 10**data1[int(idx-200000),3]/10**data1[int(idx-200000),1]
        else:
            x0 = 10**data2[int(idx-100000),3]/10**data2[int(idx-100000),1]
        ## if iii == 0:
        ##     ax.plot(xx[0], yy[0], '*', ms=7, color='k', alpha=0.3, label='At epoch of SN Ia')
        ## else:
        ##     ax.plot(xx[0], yy[0], '*', ms=7, color='k', alpha=0.3)
        out.append([xx[0], yy[0], idx])
        ##ax.annotate('%s' %idx, (xx[0], yy[0]))
    out = array(out)
    ax.set_xlabel(r'$\log$ (sSFR)')
    ax.set_ylabel(r'$\log$ (sSNR)')

    p_n = (1.19e-7, 0.586, 1.01e-11, 1.04e-9)
    p_e = (2.20e-7, 0.084, 0.55e-11, 0.41e-9)
    pcov = diag(0.2*array(p_e)**2.)
    xx = arange(-13, -7.5, 0.5)
    ax.plot(xx, log10(peices(10**xx, *p_n)), '--', color='#4BADFF',
            label='Andersen & Hjorth (2018)\n rate model')
    ps = np.random.multivariate_normal(p_n, pcov, 10000)
    ysample = asarray([log10(peices(10**xx, *pi)) for pi in ps])
    idx = where(~isnan(ysample[:,0]))
    ysample=ysample[idx]
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ## ax.fill_between(xx, upper, lower, color='red', alpha=0.2)




    ## t**(-1)
    popt=(-4.77238208e-04,  9.90040272e-01, -3.01111121e+00)
    pcov=array([[5.09893718e-08, 1.02627168e-06, 5.12290288e-06],
                [1.02627168e-06, 2.07089528e-05, 1.03610770e-04],
                [5.12290288e-06, 1.03610770e-04, 5.19460796e-04]])


    perr = sqrt(diag(pcov))
    ax.plot(xx, line_fn(xx, *popt), '--', color='#F0CBC8')#, label=r'$\beta=-1.0$')
    ps = np.random.multivariate_normal(popt, pcov, 10000)
    ysample = asarray([line_fn(xx, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)



    ## t**(-1.1)
    popt=( -1.94866195e-04,   9.95757343e-01,  -2.95934114e+00)
    pcov=array([[  8.37865874e-08,   1.68638241e-06,   8.41799528e-06],
                [  1.68638241e-06,   3.40291314e-05,   1.70253711e-04],
                [  8.41799528e-06,   1.70253711e-04,   8.53578632e-04]])
    perr = sqrt(diag(pcov))
    ax.plot(xx, line_fn(xx, *popt), '-', color='#F0CBC8', label=r'$\beta=-1.1^{+0.1}_{-0.3}$')
    ps = np.random.multivariate_normal(popt, pcov, 10000)
    ysample = asarray([line_fn(xx, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)

    ## t**(-1.4)
    popt = (-5.31848823e-05,  9.98648036e-01, -2.89739937e+00)
    pcov = array([[1.48390185e-07, 2.98669972e-06, 1.49089903e-05],
                  [2.98669972e-06, 6.02686101e-05, 3.01537222e-04],
                  [1.49089903e-05, 3.01537222e-04, 1.51179138e-03]])
    perr = sqrt(diag(pcov))
    ax.plot(xx, line_fn(xx, *popt), '--', color='#F0CBC8')#, label=r'$\beta=-1.4$')
    ps = np.random.multivariate_normal(popt, pcov, 10000)
    ysample = asarray([line_fn(xx, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)


    ## exponetial model
    popt = (0.12212169, 3.08357917, 5.71209495)
    pcov = array([[9.41823879e-05, 1.89562431e-03, 9.46250054e-03],
                  [1.89562431e-03, 3.82514562e-02, 1.91379154e-01],
                  [9.46250054e-03, 1.91379154e-01, 9.59494348e-01]])

    perr = sqrt(diag(pcov))
    ax.plot(xx, line_fn(xx, *popt), 'b-', label=r'Exponential model')
    ps = np.random.multivariate_normal(popt, pcov, 10000)
    ysample = asarray([line_fn(xx, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(xx, upper, lower, color='green', alpha=0.2)
    lower = percentile(ysample, 2.5, axis=0)
    upper = percentile(ysample, 97.5, axis=0)
    ax.fill_between(xx, upper, lower, color='green', alpha=0.2)

    p0 = (1.,1.,0.0)
    popt, pcov = curve_fit(line_fn, out[:,0], out[:,1], p0=p0)
    perr = sqrt(diag(pcov))
    #ax.plot(xx, line_fn(xx, *popt), 'g-', label=r'$\Phi(t)\propto$G(2.5,1..0)')
    ps = np.random.multivariate_normal(popt, pcov, 10000)
    ysample = asarray([line_fn(xx, *pi) for pi in ps])
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    #ax.fill_between(xx, upper, lower, color='green', alpha=0.2)
    print(popt,pcov)

    
    ## plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    ## cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    ## cbar = plt.colorbar(j1, cax=cax)
    ## cbar.ax.set_yticklabels(arange(0, 15, 1))
    ## cbar.ax.set_ylabel('Lookback Time (Gyr)')

    ax2 = ax.twinx()
    data = loadtxt('mannucci05.txt')
    ax2.errorbar(data[1:,0], log10(data[1:,1]),
                 yerr=log10((data[1:,1]-data[1:,2])/data[1:,1]),
                 fmt='o', ms=10, mfc='white', label='Mannucci et al. (2005)')
    ax2.errorbar(data[0,0], log10(data[0,1]),
                 yerr=0.5,xerr=data[0,0]+11,
                 uplims=0.5,
                 fmt='o', ms=10, mfc='white', color='#2077B4')



    data = loadtxt('sullivan06.txt')
    ax2.errorbar(data[1:,0], log10(data[1:,1]),
                 yerr=-log10((data[1:,1]-data[1:,2])/data[1:,1]),
                 fmt='o', ms=10, mfc='white', label='Sullivan et al. (2006)')
    ax2.errorbar(data[0,0], log10(data[0,1]),
                 yerr=0.5,xerr=data[0,0]+11,
                 uplims=0.5,
                 fmt='o', ms=10, mfc='white',color='#FF7F0E')

    data = loadtxt('smith12.txt')
    lolims = zeros((len(data)),)
    lolims[0]=1.0
    ax2.errorbar(data[:,0], log10(data[:,1]),
                 yerr=[-log10((data[:,1]-data[:,2])/data[:,1]),log10((data[:,1]-data[:,2])/data[:,1])],
                 uplims=lolims,
                 fmt='o', ms=10, mfc='white', label='Smith et al. (2012)')
    ax2.errorbar(data[0,0], log10(data[0,1]),
                 yerr=0.5,xerr=data[0,0]+11,
                 uplims=0.5,
                 fmt='o', ms=10, mfc='white', color='#2CA02C')
    

    ax.set_xticks(list(arange(-12.5, -7,1)))
    ax.set_yticks(list(arange(-14, -10, 1)))
    ax2.set_yticks([])
    ax.set_xlim(-13, -8)
    ax.set_ylim(-14.25, -11.25)
    ax2.set_ylim(-14.25, -11.25)

    ax.legend(loc=2,frameon=False)
    ax2.legend(loc=4,frameon=False)

    tx = arange(-13,-10, 1) 
    #ax.fill_between(tx, log10(p_n[0]*p_n[2]**p_n[1]), -15, color='0.75', alpha=0.3)
    tx = arange(-9,-7, 1) 
    #ax.fill_between(tx, log10(p_n[0]*p_n[3]**p_n[1]), -11, color='0.75', alpha=0.3)


    ## ax2.annotate('Passive Gal. & UDGs', xy=(-12.4,-13.8), xycoords='data',
    ##             rotation=0, size=10, zorder=100,
    ##             bbox=dict(boxstyle="round4", fc="w", ec='w', alpha=0.7)
    ##             )
    
    ## ax.annotate('Green\n Peas', (-8.75,-11.75), xycoords='data',
    ##             rotation=0, size=10, zorder=100,
    ##             bbox=dict(boxstyle="round4", fc="w", ec='w', alpha=0.7)
    ##             )
                

    savefig('figure_ssfr.png')
    
            
        
        
        
        
            
