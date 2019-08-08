#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from scipy.stats import ks_2samp as ks

rcParams['figure.figsize']=20, 10
rcParams['font.size']=24.0

file1 = '../ALLSFH_new_z/Cami_GOODS-N_zbest.dat'
file2 = '../ALLSFH_new_z/Cami_GOODS-S_zbest.dat'

data1 = loadtxt(file1)
data1[:,0] = data1[:,0]+100000
data2 = loadtxt(file2)
data2[:,0] = data2[:,0]+200000
data = concatenate((data1,data2), axis=0)

idx = loadtxt('../DTD_fits/host_idxs.txt')

ax = subplot(121) #Mass distributions
tmp = data[:,1][where((data[:,1]>-99)&(data[:,1]!=8.0))]
## line = u.binmode(tmp, bins=200)[0]
N1, bins1, j1 = ax.hist(tmp, bins=50, color='green', histtype='step', cumulative=True, normed=True)

mask = in1d(data[:,0], idx)
data = data[mask]

#tmp = data[:,1][where(data[:,1]>-99)]
tmp = data[:,1][where((data[:,1]>-99)&(data[:,1]!=8.0))]
## line = u.binmode(tmp,bins=6)[0]
N2, bins2, j2 = ax.hist(tmp, color='orange', cumulative=True, normed=True, alpha=0.3)

dval, pval=ks(N1, N2)
print(dval, pval)

conf = 0.95
ca = sqrt(-1/2.*log(1-conf))
dnp = ca*sqrt((len(bins1)+len(bins2)-2)/((len(bins1)-1)*(len(bins2)-1)))
print(dnp)
if dval <= dnp:
    print('Null hypothesis cannot be rejected at %d level' %(int(conf*100.)))
else:
    print('Distributions are different at %d level' %(int(conf*100.)))


ax2 = subplot(122) #SFR
tmp = data[:,3][where((data[:,3]>-99))]#&(data[:,1]!=8.0))]
## line = u.binmode(tmp, bins=200)[0]
N1, bins1, j1 = ax2.hist(tmp, bins=200, color='green', histtype='step', cumulative=True, normed=True, label='All Catalog Galaxies')

mask = in1d(data[:,0], idx)
data = data[mask]

tmp = data[:,3][where(data[:,3]>-99)]
#tmp = data[:,1][where((data[:,1]>-99)&(data[:,1]!=8.0))]
## line = u.binmode(tmp,bins=6)[0]
N2, bins2, j2 = ax2.hist(tmp, color='orange', cumulative=True, normed=True, alpha=0.3, label='SN Ia Hosts')

dval, pval=ks(N1, N2)
print(dval, pval)

conf = 0.995
ca = sqrt(-1/2.*log(1-conf))
dnp = ca*sqrt((len(bins1)+len(bins2)-2)/((len(bins1)-1)*(len(bins2)-1)))
print(dnp)
if dval <= dnp:
    print('Null hypothesis cannot be rejected at %d level' %(int(conf*100.)))
else:
    print('Distributions are different at %d level' %(int(conf*100.)))

ax.set_xlabel(r'Log(M/M$_\odot$)')
ax2.set_xlabel('Log(SFR)')
ax.set_ylabel('CDF')

u.adjust_spines(ax, (['left','bottom']))
u.adjust_spines(ax2, (['bottom']))
ax.set_xlim(8,11.6)
ax2.set_xlim(-2,2.5)
ax2.legend(loc=2, frameon=False)
tight_layout()
savefig('figure_cdfs.png')
