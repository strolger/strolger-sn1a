#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import pickle
rcParams['font.size']=12.0
#rcParams['figure.figsize']=8,8

def redchisqg(ydata,ymod,sd, deg=2):  
      # Chi-square statistic  
      chisq=sum( ((ydata-ymod)/sd)**2 )  
             
      # Number of degrees of freedom assuming 2 free parameters  
      nu=ydata.size-1-deg  
        
      return chisq/nu       


scale = (0.0210)*0.062*(0.7)**2
ax = subplot(111)
ax2 = ax.twinx()

p0 = [-1257.93, 59.32, 248.88] ## from MCMC


time = arange(0.05,13.5,0.1)
dtd = rz.dtdfunc(time,*p0)

def func(x,*p):
    aa=p
    return(aa*rz.dtdfunc(x,*p0)*scale)
p_t = [1.0]


color = '#4BADFF'
## pwrl = (-1.0,1.0)
## ax.plot(time,rz.powerdtd(time, *pwrl, normed=False)*(1.3e-3), 'b--', label=r'$t^{%.1f}$'%(pwrl[0]))
pwrl = (-1.10,1.0)
perr = (0.09, 0.0)
cov = diag(array(perr)**2)
ps = np.random.multivariate_normal(pwrl, cov, 100)
ysample = asarray([rz.powerdtd(time, *pi, normed=False)*(1.6e-3) for pi in ps])
lower = percentile(ysample, 15.9, axis=0)
upper = percentile(ysample, 84.1, axis=0)
ax.fill_between(time, upper, lower, color=color, alpha=0.2)
ax.plot(time,rz.powerdtd(time, *pwrl, normed=False)*(1.6e-3), '--', color=color, label=r'Field ($\beta={%.2f}^{+0.08}_{-0.07}$)'%(pwrl[0]))



color = '#F0CBC8'
pwrl = (-1.39,1.0)
perr = (0.32, 0.0)
cov = diag(array(perr)**2)
ps = np.random.multivariate_normal(pwrl, cov, 100)
ysample = asarray([rz.powerdtd(time, *pi, normed=False)*(5.4e-3) for pi in ps])
lower = percentile(ysample, 15.9, axis=0)
perr = (0.05, 0.0)
cov = diag(array(perr)**2)
ps = np.random.multivariate_normal(pwrl, cov, 100)
ysample = asarray([rz.powerdtd(time, *pi, normed=False)*(5.4e-3) for pi in ps])
upper = percentile(ysample, 84.1, axis=0)
ax.fill_between(time, upper, lower, color=color, alpha=0.2)
ax.plot(time,rz.powerdtd(time, *pwrl, normed=False)*(5.4e-3), '--', color=color, label=r'Clusters ($\beta={%.1f}^{+0.32}_{-0.05}$)'%(pwrl[0]))


data = loadtxt('maoz_table2.txt')
msc = 1.0e-4*0.7**2
idx = where(data[:,6]==1)
ax2.errorbar(data[idx][:,0], data[idx][:,3]*msc,
            xerr=[-1*data[idx][:,1],data[idx][:,2]],
            yerr=[-1*data[idx][:,4]*msc,data[idx][:,5]*msc],
            fmt='o', color='k', mfc='1.0', mec='red', label='Cluster Rates')
pwrl = (-1.39,1.0)
tmp=rz.powerdtd(data[idx][:,0], *pwrl, normed=False)*(5.4e-3)
chi2 = redchisqg(data[idx][:,3]*msc, tmp, sd=data[idx][:,5]*msc)
print(chi2, len(tmp)-1)


idx = where(data[:,6]==0)

ax2.errorbar(data[idx][1:,0], data[idx][1:,3]*msc,
            xerr=[-1*data[idx][1:,1],data[idx][1:,2]],
            yerr=[-1*data[idx][1:,4]*msc,data[idx][1:,5]*msc],
            fmt='o', color='k', mfc='1.0', mec='blue', label='Field Rates')
ax2.errorbar(data[idx][0,0], data[idx][0,3]*msc,
             xerr=data[idx][0,2],
             yerr=0.05,#data[idx][0,5]*msc,#,data[idx][:,5]*msc],
             lolims=data[idx][0,5]*msc,
             fmt='o', color='k', mfc='1.0', marker='^', ms=10, mec='blue', label='LMC & SMC limit')
pwrl = (-1.10,1.0)
tmp=rz.powerdtd(data[idx][1:,0], *pwrl, normed=False)*(1.6e-3)
chi2 = redchisqg(data[idx][1:,3]*msc, tmp, sd=data[idx][1:,5]*msc, deg=0)
print(chi2, len(tmp)-1)






#popt, pcov = curve_fit(func, data[:,0], data[:,3]*msc, p0=p_t, sigma=data[:,4]*msc)
#print(popt)
popt = [12.921]
ax.plot(time,dtd*scale*popt,'b-', lw=2, label= 'Exponential model')#label='Norm = %1.1f' %(simps(dtd,x=time)))
tmp=rz.dtdfunc(data[:,0], *p0)*scale*popt
chi2 = redchisqg(data[:,3]*msc, tmp, sd=data[:,5]*msc, deg=0)
print(chi2,len(tmp)-1)

files = glob.glob('mc_sfd_*.pkl')
tot = 0
for file in files:
    samples = pickle.load(open(file,'rb'))
    samples = samples[:,50:,:].reshape((-1,5))
    print('adding %d samples from %s... '%(len(samples), file))
    try:
        temp = concatenate((temp, samples), axis=0)
    except:
        temp = samples

print('%d total samples' %len(temp))
samples = temp
ndraws = 100 #21 is actually pretty good
draws  = samples[np.random.randint(0, len(samples), ndraws),:]
ysample = asarray([rz.dtdfunc(time, *pi[1:-1])*exp(pi[0])*(0.0210)*(0.7)**2*popt for pi in draws])


lower = percentile(ysample, 15.9, axis=0)
upper = percentile(ysample, 84.1, axis=0)
ax.fill_between(time, upper, lower, color='green', alpha=0.2)
lower = percentile(ysample, 2.5, axis=0)
upper = percentile(ysample, 97.5, axis=0)
ax.fill_between(time, upper, lower, color='green', alpha=0.2)


ax.set_ylabel(r'SN Ia yr$^{-1}$ ($10^{10}$ M$_{\odot}$)$^{-1}$')
ax.set_xlabel('Delay Time (Gyr)')
ax.set_xlim(0.03, 14)
ax.set_ylim(5e-5,0.1)
ax2.set_ylim(5e-5,0.1)
ax.legend(loc=1,numpoints=1,frameon=False)
ax2.legend(loc=3,numpoints=1,frameon=False)
ax.set_xscale('log')
ax.set_yscale('log')
ax2.set_yscale('log')
ax2.set_yticks([])
savefig('figure_loglog_dtd.png')
##show()
