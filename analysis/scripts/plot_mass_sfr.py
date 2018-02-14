#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u





if __name__=='__main__':
    

    file = sys.argv[1]
    f = open(file)
    lines = f.readlines()
    f.close()
    data = []
    for line in lines:
        if line.startswith('#'): continue
        c = line.split()
        data.append(c[5:8])
    data = array(data)

    ax1=subplot(111)
    ax1.plot(data[:,0],data[:,1],'k.')

    ax1.set_xlabel('LogM_50')
    ax1.set_ylabel('LogSFR_50')

    show()
    
                 
