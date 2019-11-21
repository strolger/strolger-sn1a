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
allmass = data[:,1][where((data[:,1]>-99)&(data[:,1]!=8.0))]
N1, bins1, j1 = ax.hist(allmass, bins=100, color='green', histtype='step', cumulative=True, normed=True)
savefile = False
if savefile:
    outfile = 'all_masses.txt'
    if os.path.isfile(outfile): os.remove(outfile)
    f = open(outfile,'w')
    for item in allmass:f.write('%s\n'%item)
    f.close()


mask = in1d(data[:,0], idx)
data_ias = data[mask]

hostmass = data_ias[:,1][where((data_ias[:,1]>-99)&(data_ias[:,1]!=8.0))]
N2, bins2, j2 = ax.hist(hostmass, bins=20, color='orange', cumulative=True, normed=True, alpha=0.3)
if savefile:
    outfile = 'host_masses.txt'
    if os.path.isfile(outfile): os.remove(outfile)
    f = open(outfile,'w')
    for item in hostmass:f.write('%s\n'%item)
    f.close()

dval, pval=ks(hostmass, allmass)
print(dval, pval)

conf = 0.99
ca = sqrt(-1/2.*log(1.-conf))
print(ca)
dnp = ca*sqrt((len(hostmass)+len(allmass)-2)/((len(hostmass)-1)*(len(allmass)-1)))
print(dnp)
if dval <= dnp:
    print('Null hypothesis cannot be rejected at %d%% level' %(int((conf)*100.)))
else:
    print('Distributions are different at %d%% level' %(int((conf)*100.)))

ax2 = subplot(122) #SFR
allsfr = data[:,3][where((data[:,3]>-99))]
N1, bins1, j1 = ax2.hist(allsfr, bins=100, color='green', histtype='step', cumulative=True, normed=True, label='All Catalog Galaxies')
if savefile:
    outfile = 'all_sfrs.txt'
    if os.path.isfile(outfile): os.remove(outfile)
    f = open(outfile,'w')
    for item in tmp:f.write('%s\n'%item)
    f.close()

hostsfr = data_ias[:,3][where(data_ias[:,3]>-99)]
N2, bins2, j2 = ax2.hist(hostsfr, bins=20, color='orange', cumulative=True, normed=True, alpha=0.3, label='SN Ia Hosts')
if savefile:
    outfile = 'host_sfrs.txt'
    if os.path.isfile(outfile): os.remove(outfile)
    f = open(outfile,'w')
    for item in tmp:f.write('%s\n'%item)
    f.close()

dval, pval=ks(hostsfr,allsfr)
print(dval, pval)

conf = 0.99
ca = sqrt(-1/2.*log(1-conf))
dnp = ca*sqrt((len(hostsfr)+len(allsfr)-2)/((len(hostsfr)-1)*(len(allsfr)-1)))
print(dnp)
if dval <= dnp:
    print('Null hypothesis cannot be rejected at %d%% level' %(int(conf*100.)))
else:
    print('Distributions are different at %d%% level' %(int(conf*100.)))

ax.set_xlabel(r'Log(M/M$_\odot$)')
ax2.set_xlabel('Log(SFR)')
ax.set_ylabel('CDF')

u.adjust_spines(ax, (['left','bottom']))
u.adjust_spines(ax2, (['bottom']))
ax.set_xlim(8,11.7)
ax2.set_xlim(-2,2.8)
ax2.legend(loc=2, frameon=False)
tight_layout()
savefig('figure_cdfs.png')
