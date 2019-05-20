#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u

if __name__=='__main__':
    
    data = loadtxt('CANDELS_GDSN_znew_avgal_radec.dat')

    ax = subplot(111)

    phot = data[:,8]
    phot[where(phot>90)]=float('nan')
    ax.hist(phot[~isnan(phot)])

    show()
    
