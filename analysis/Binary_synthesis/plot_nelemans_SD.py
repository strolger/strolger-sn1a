#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from scipy.optimize import curve_fit


if __name__=='__main__':


    fig= plt.figure()
    grid = plt.GridSpec(1,3,wspace=0.2)

    data = loadtxt('nelemans_table.txt')

    steps = False

    ax = fig.add_subplot(grid[0,:-1])
    ax2 = fig.add_subplot(grid[0,2])
    tt = arange(0,10000,100)
    p0 = [
        [10, 600, 220],
        [110, 1000, 2],
        [350, 1200, 20],
        [6000, 6000, -2],
        ]
    scales = [
        6.,
        2.,
        6.,
        8e-2,
        ]
        
    ## for ii in range(shape(data)[1]-2):
    ##     ax.errorbar(data[:,0], data[:,ii+2]*1e-13, xerr=data[:,1], fmt='o')
    ##     ax.plot(tt, rz.dtdfunc(tt, *p0[ii])*scales[ii])

    ii=0
    if steps:
        ax.plot (data[:,0]+data[:,1], data[:,2]*1e-3, drawstyle='steps', color=u.my_color(0), alpha=0.3, label = 'Mennekens et al. 2010')
    else:
        ax.errorbar(data[:,0], data[:,2]*1e-3, xerr=data[:,1], fmt='o',color=u.my_color(0), label='Mennekens et al. 2010')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(0))
    ii=1
    if steps:
        ax.plot (data[:,0]+data[:,1], data[:,3]*1e-3, drawstyle='steps', color=u.my_color(2), alpha=0.3, label = 'Yugelson 2010')
    else:
        ax.errorbar(data[:,0], data[:,3]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(2), label='Yungelson 2010')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(2))
    ii=2
    if steps:
        ax.plot (data[:,0]+data[:,1], data[:,5]*1e-3, drawstyle='steps', color=u.my_color(5), alpha=0.3, label = 'Wang et al. 2010')
    else:
        ax.errorbar(data[:,0], data[:,5]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(5), label='Wang et al. 2010')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(5))
    ii=3
    if steps:
        ax.plot (data[:,0]+data[:,1], data[:,6]*1e-3, drawstyle='steps', color=u.my_color(6), alpha=0.3, label = 'Ruiter et al. 2009')
    else:
        ax.errorbar(data[:,0], data[:,6]*1e-3, xerr=data[:,1], fmt='o', color=u.my_color(6), label='Ruiter et al. 2009')
    ax.plot(tt,rz.dtdfunc(tt,*p0[ii])*scales[ii], color=u.my_color(6))


    data = loadtxt('nelemans_table2.txt')
    if steps:
        ax2.plot (data[:,0]+data[:,1], data[:,2]*1e-3, drawstyle='steps', color=u.my_color(2), label = 'Yungelson 2010')
    else:
        ax2.errorbar(data[:,0], data[:,2]*1e-3, xerr=data[:,1], fmt='o',color=u.my_color(2), label='Yungelson 2010')

    p1=[[-650, 2200, 1100]]
    ax2.plot(tt, 5*rz.dtdfunc(tt,*p1[0]), '-', color=u.my_color(2))#, label='best model fit')
    pa = (-1.0, 1.0)
    yy = rz.powerdtd(tt, *pa)
    ax2.plot(tt, yy*2.0,'--', color=u.my_color(2), label=r'$t^{-1}$')


    col_labels = [r'$\xi$', r'$\omega$', r'$\alpha$']
    table_vals = p0
    the_table = ax.table(cellText=table_vals,
                         colWidths = [0.1]*3,
                         #rowLabels=row_labels,
                         colLabels=col_labels,
                         cellLoc='center',
                         loc='best')
    ax.text(12,3.4,'Table Title',size=9)
    for key, cell in the_table.get_celld().items():
        cell.set_linewidth(0)

    table_vals = p1
    the_table = ax2.table(cellText=table_vals,
                          colWidths = [0.26]*3,
                          #rowLabels=row_labels,
                          colLabels=col_labels,
                          cellLoc='center',
                          loc='best')
    ax2.text(12,3.4,'Table Title',size=9)
    for key, cell in the_table.get_celld().items():
        cell.set_linewidth(0)


    #ax2.set_xscale('log')
    #ax2.set_yscale('log')

    ax2.set_xlim(100,10000)
    ## ax2.set_ylim(1e-5, 0.005)
    ax2.set_ylim(4e-7, 2e-2)
    ax2.set_xlabel('Delay time (Myr)')


    ax2.set_yticks([])
    #ax2.legend(loc=3, fontsize=8,frameon=False)
    ax2.set_title('DD Models')

    
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_ylim(4e-7,1.25e-2)
    #ax.legend(loc=3, frameon=True)

    ax.set_xlabel('Delay time (Myr)')
    ax.set_ylabel(r'SN Ia yr$^{-1}$ ($10^{10}$ M$_{\odot}$)$^{-1}$')
    ax.set_title('SD Models')
    savefig('figure_sd_dd_fits.pdf')
    
    
    
