#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u



file1 = '../ALLSFH_new_z/Cami_GOODS-N_zbest.dat'
file2 = '../ALLSFH_new_z/Cami_GOODS-S_zbest.dat'

data1 = loadtxt(file1)
data1[:,0] = data1[:,0]+100000
data2 = loadtxt(file2)
data2[:,0] = data2[:,0]+200000
data = concatenate((data1,data2), axis=0)

idx = loadtxt('host_idxs.txt')

ax = subplot(111)
tmp = data[:,1][where((data[:,1]>-99)&(data[:,1]!=8.0))]
line = u.binmode(tmp, bins=200)[0]
N, bins, j1 = ax.hist(tmp, bins=200, color='green')
ax.axvline(line, color='b')
print(line, average(tmp), std(tmp))

ax2 = ax.twinx()

mask = in1d(data[:,0], idx)
data = data[mask]

tmp = data[:,1][where(data[:,1]>-99)]
## tmp = data[:,1][where((data[:,1]>-99)&(data[:,1]!=8.0))]
line = u.binmode(tmp,bins=6)[0]
N, bins, j1 = ax2.hist(tmp, bins=6, color='orange', alpha=0.3)
ax2.axvline(line, color='r')
print(line, average(tmp), std(tmp))

ax2.set_ylim(0,100)

savefig('figure_avemasses.png')
