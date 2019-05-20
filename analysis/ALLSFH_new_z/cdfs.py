#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u

if __name__=='__main__':
    

    ax=subplot(111)

    d1 = loadtxt('Cami_GOODS-N_zbest.dat')
    d2 = loadtxt('Cami_GOODS-S_zbest.dat')

    hosts = loadtxt('../DTD_fits/host_idxs.txt')

    data = concatenate((d1, d2), axis=0)
    blogmass = data[:,1]
    blogmass[where(blogmass==-99.0)]=float('nan')
    ax.hist(blogmass[~isnan(blogmass)])

    

    show()
    

    
