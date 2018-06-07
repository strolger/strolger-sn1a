#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz


if __name__=='__main__':
    

    p0 = (3.5, 2.5, 2.5)
    tt = arange(-10.0, 13.6, 0.1)
    ax1 = subplot(131)
    for mm in arange(-2.5, 12.5, 3):
        p = (mm,p0[1],p0[2])
        ax1.plot(tt, rz.dtdfunc(tt,*p, norm=False))
    ax1.set_xlim(0,10)
    ax2 = subplot(132)
    for ww in [1.0, 2.0, 2.5, 5.0, 7.5]:
        p = (p0[0],ww,p0[2])
        ax2.plot(tt, rz.dtdfunc(tt,*p, norm=False))
    ax2.set_xlim(0,10)
    ax3 = subplot(133)
    for kk in arange(-5.0, 7.5, 2.5):
        p = (p0[0],p0[1],kk)
        ax3.plot(tt, rz.dtdfunc(tt,*p, norm=False))
    ax3.set_xlim(0,10)


    u.adjust_spines(ax1,['left','bottom'])
    u.adjust_spines(ax2,['bottom'])
    u.adjust_spines(ax3,['bottom'])
    
    ax1.set_ylabel(r'$\Phi(\tau)$')
    ax2.set_xlabel('Time (Gyr)')
    tight_layout()
    savefig('dtd_families.png')
    
    
