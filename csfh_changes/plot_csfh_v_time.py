#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import cosmotools as ct



if __name__=='__main__':


    p0 = [0.013, 2.6, 3.2, 6.1]
    redshifts = arange(0, 6.1, 0.1)
    lbt = array([ct.cosmotime(x) for x in redshifts])
    tt = 13.6 - lbt
    csfh = rz.csfh(redshifts, *p0)
    ax = subplot(111)

    ax.plot(tt, csfh, 'r-')

    tt = arange(0, 13.6, 0.05)
    sfh_t = rz.sfr_behroozi_12(tt)
    ax.plot(tt,sfh_t, 'b-')

    show()
    
