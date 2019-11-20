#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u

    
all_mass = loadtxt('all_masses.txt')
host_mass = loadtxt('host_masses.txt')

ax1 = subplot(121)

N1, bins, j1 = ax1.hist(host_mass, normed=False,zorder=2, color='C1', label='Host Galaxies')
N2, bins, j2 = ax1.hist(all_mass, bins=bins, normed=False,zorder=1, color='C0', label='All Galaxies')

ax1.set_yscale('log')
ax1.set_xlabel('Log(M$_\odot$)')
ax1.set_ylabel('Number')

ax1.legend(loc=1)#, frameon=False)


ax2 = subplot(122)

all_sfr = loadtxt('all_sfrs.txt')
host_sfr = loadtxt('host_sfrs.txt')

N3, bins, j3 = ax2.hist(host_sfr, normed=False, zorder=2, color='C1')
N4, bins, j4 = ax2.hist(all_sfr, bins=bins, normed=False, zorder=1, color='C0')

ax2.set_yscale('log')
ax2.set_xlabel('Log(SFR)')
#ax1.set_ylabel('Number')

savefig('figure_host_distribution.png')



