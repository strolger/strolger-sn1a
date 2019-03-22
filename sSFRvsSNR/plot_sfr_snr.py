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


def highResPoints(x,y,factor=10):
    '''
    Take points listed in two vectors and return them at a higher
    resultion. Create at least factor*len(x) new points that include the
    original points and those spaced in between.

    Returns new x and y arrays as a tuple (x,y).
    '''

    # r is the distance spanned between pairs of points
    r = [0]
    for i in range(1,len(x)):
        dx = x[i]-x[i-1]
        dy = y[i]-y[i-1]
        r.append(np.sqrt(dx*dx+dy*dy))
    r = np.array(r)

    # rtot is a cumulative sum of r, it's used to save time
    rtot = []
    for i in range(len(r)):
        rtot.append(r[0:i].sum())
    rtot.append(r.sum())

    dr = rtot[-1]/(NPOINTS*RESFACT-1)
    xmod=[x[0]]
    ymod=[y[0]]
    rPos = 0 # current point on walk along data
    rcount = 1 
    while rPos < r.sum():
        x1,x2 = x[rcount-1],x[rcount]
        y1,y2 = y[rcount-1],y[rcount]
        dpos = rPos-rtot[rcount] 
        theta = np.arctan2((x2-x1),(y2-y1))
        rx = np.sin(theta)*dpos+x1
        ry = np.cos(theta)*dpos+y1
        xmod.append(rx)
        ymod.append(ry)
        rPos+=dr
        while rPos > rtot[rcount+1]:
            rPos = rtot[rcount+1]
            rcount+=1
            if rcount>rtot[-1]:
                break

    return xmod,ymod


## #CONSTANTS
## NPOINTS = 10
## COLOR='blue'
## RESFACT=10
## MAP='winter' # choose carefully, or color transitions will not appear smoooth

## # create random data
## np.random.seed(101)
## x = np.random.rand(NPOINTS)
## y = np.random.rand(NPOINTS)

## fig = plt.figure()
## ax1 = fig.add_subplot(221) # regular resolution color map
## ax2 = fig.add_subplot(222) # regular resolution alpha
## ax3 = fig.add_subplot(223) # high resolution color map
## ax4 = fig.add_subplot(224) # high resolution alpha

## # Choose a color map, loop through the colors, and assign them to the color 
## # cycle. You need NPOINTS-1 colors, because you'll plot that many lines 
## # between pairs. In other words, your line is not cyclic, so there's 
## # no line from end to beginning
## cm = plt.get_cmap(MAP)
## ax1.set_color_cycle([cm(1.*i/(NPOINTS-1)) for i in range(NPOINTS-1)])
## for i in range(NPOINTS-1):
##     ax1.plot(x[i:i+2],y[i:i+2])


## ax1.text(.05,1.05,'Reg. Res - Color Map')
## ax1.set_ylim(0,1.2)

## # same approach, but fixed color and 
## # alpha is scale from 0 to 1 in NPOINTS steps
## for i in range(NPOINTS-1):
##     ax2.plot(x[i:i+2],y[i:i+2],alpha=float(i)/(NPOINTS-1),color=COLOR)

## ax2.text(.05,1.05,'Reg. Res - alpha')
## ax2.set_ylim(0,1.2)

## # get higher resolution data
## xHiRes,yHiRes = highResPoints(x,y,RESFACT)
## npointsHiRes = len(xHiRes)

## cm = plt.get_cmap(MAP)

## ax3.set_color_cycle([cm(1.*i/(npointsHiRes-1)) 
##                      for i in range(npointsHiRes-1)])


## for i in range(npointsHiRes-1):
##     ax3.plot(xHiRes[i:i+2],yHiRes[i:i+2])

## ax3.text(.05,1.05,'Hi Res - Color Map')
## ax3.set_ylim(0,1.2)

## for i in range(npointsHiRes-1):
##     ax4.plot(xHiRes[i:i+2],yHiRes[i:i+2],
##              alpha=float(i)/(npointsHiRes-1),
##              color=COLOR)
## ax4.text(.05,1.05,'High Res - alpha')
## ax4.set_ylim(0,1.2)



## fig.savefig('gradColorLine.png')
## plt.show()




if __name__=='__main__':

    if not os.path.isfile('ssSFRN.pkl'):
        idxs = loadtxt('../analysis/scripts/host_idxs.txt')
        sfhs = {}
        for idx in idxs:
            print(idx)
            if idx > 200000:
                file = '../analysis/ALLSFH_new_z/gsznnpas/%05d.dat'%(idx-200000)
            else:
                file = '../analysis/ALLSFH_new_z/gnznnpas/%05d.dat'%(idx-100000)
            if os.path.isfile(file):
                data = loadtxt(file)
            else:
                print('%s not found'%file)
            ndata = get_data(data, p0=(-1258, 59, 248), frac_ia = 0.06, plotit=True)
            sfhs[idx]=ndata
            pickle.dump(sfhs,open('ssSFRN.pkl','wb'))
    else:
        sfhs = pickle.load(open('ssSFRN.pkl','rb'))
        

    data1 = loadtxt('../analysis/ALLSFH_new_z/Cami_GOODS-S_zbest.dat')
    data2 = loadtxt('../analysis/ALLSFH_new_z/Cami_GOODS-N_zbest.dat')
    

    ax = subplot(111)
    k1 = list(sfhs.keys())[0]
    ii = max(where(sfhs[k1][:,0]<0.05)[0])
    for idx in sfhs.keys():
        xx = log10(sfhs[idx][ii:,1])
        yy = log10(sfhs[idx][ii:,2])
        y0 = sfhs[idx][0,2]
        j1 = ax.scatter(xx, yy, c=range(len(xx)), cmap='plasma', marker='_', s=1)

        ## cm = plt.get_cmap('bone')
        ## ax.set_color_cycle([cm(1.*i/(len(xx)-1)) for i in range(len(xx)-1)])
        ## for i in range(len(xx)-1):
        ##     ax.plot(xx[i:i+2], yy[i:i+2])


        if idx > 200000:
            x0 = 10**data1[int(idx-200000),3]/10**data1[int(idx-200000),1]
        else:
            x0 = 10**data2[int(idx-100000),3]/10**data2[int(idx-100000),1]
        ## ax.plot(log10(x0), log10(y0), '*', color='k')
        ax.plot(xx[0], yy[0], '*', ms=7, color='k')
        ## ax.annotate('%s' %idx, (xx[0], yy[0]))
    ax.set_xlabel(r'$\log$ (sSFR)')
    ax.set_ylabel(r'$\log$ (sSNR)')


    p_n = (1.19e-7, 0.586, 1.01e-11, 1.04e-9)
    p_e = (2.20e-7, 0.084, 0.55e-11, 0.41e-9)
    pcov = diag(array(p_e)**2.)
    xx = arange(-13, -7.5, 0.5)
    ax.plot(xx, log10(peices(10**xx, *p_n)), 'k-')
    ps = np.random.multivariate_normal(p_n, pcov, 10000)
    ysample = asarray([log10(peices(10**xx, *pi)) for pi in ps])
    idx = where(~isnan(ysample[:,0]))
    ysample=ysample[idx]
    lower = percentile(ysample, 15.9, axis=0)
    upper = percentile(ysample, 84.1, axis=0)
    ax.fill_between(xx, upper, lower, color='red', alpha=0.2)
  



    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    
    
    cbar = plt.colorbar(j1, cax=cax)
    cbar.ax.set_yticklabels(arange(0, 15, 1))
    cbar.ax.set_ylabel('Lookback Time (Gyr)')
    savefig('figure_ssfr.png')
    
            
        
        
        
        
            
