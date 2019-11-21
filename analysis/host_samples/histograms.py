#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u

rcParams['figure.figsize']=20, 10
rcParams['font.size']=24.0

all_mass = loadtxt('all_masses.txt')
host_mass = loadtxt('host_masses.txt')

ax1 = subplot(121)

N2, bins, j2 = ax1.hist(all_mass, normed=False,zorder=1, color='C0', label='All Catalog Galaxies')
N1, bins, j1 = ax1.hist(host_mass, bins=bins, normed=False, zorder=2, color='C1', label='SN Ia Hosts')



ax2 = subplot(122)

all_sfr = loadtxt('all_sfrs.txt')
host_sfr = loadtxt('host_sfrs.txt')

N3, bins, j3 = ax2.hist(host_sfr, normed=False, zorder=2, color='C1')
N4, bins, j4 = ax2.hist(all_sfr, bins=bins, normed=False, zorder=1, color='C0')


ax1.set_yscale('log')
ax1.set_xlabel('Log(M/M$_\odot$)')
ax1.set_ylabel('Number')

ax1.legend(loc=1, frameon=False)

ax2.set_yscale('log')
ax2.set_xlabel('Log(SFR)')
#ax1.set_ylabel('Number')

u.adjust_spines(ax1, (['left','bottom']))
u.adjust_spines(ax2, (['right','bottom']))
ax1.set_xlim(8,11.7)
ax2.set_xlim(-2,2.8)
tight_layout()
savefig('figure_host_distribution.png')



