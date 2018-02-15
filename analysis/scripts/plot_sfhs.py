#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from scipy.integrate import simps
from stsci.convolve import boxcar
rcParams['figure.figsize']=15.3, 18.75
rcParams['font.size']=18.0

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


if __name__=='__main__':



    sfhs = glob.glob('../SN_SFH/*/*.dat')
    lbu = 9.5
    lbl = 0.05
    bxsmooth=3
    
    fig = plt.figure()
    fig.suptitle = ('test')
    fig.text(0.5,0.04,'Lookback Time (Gyr)',ha='center',va='center')
    fig.text(0.06,0.5,'Normalized Star Formation Rate',ha='center',va='center',rotation='vertical')
    

    
    ax1=subplot(221)
    composite = 0
    for sfh in sfhs:
        data = loadtxt(sfh)[1:]
        ii = where((data[:,0]>lbl)&(data[:,0]<=lbu)) # cuts the data range
        data = data[ii]
        data[:,1] = boxcar(data[:,1],(bxsmooth,))
        composite+=data[:,1]/simps(data[:,1])
        ax1.plot(data[:,0], data[:,1]/simps(data[:,1])/len(sfhs)*30, '-', color='0.65',alpha=0.3)
    ax1.plot(data[:,0],composite/len(sfhs),'k-', label='%d SNe Hosts' %len(sfhs))

    codex = []
    f = open('../SN_SFH/goodsn_match_idcand.dat')
    lines=f.readlines()
    f.close()
    for line in lines:
        if line.startswith('#'): continue
        c = line.split()
        codex.append(c[:2])
    f = open('../SN_SFH/goodss_match_idcand.dat')
    lines=f.readlines()
    f.close()
    for line in lines:
        if line.startswith('#'): continue
        c = line.split()
        codex.append(c[:2])
    codex=array(codex)

    f = open('goods_1as.csv')
    lines = f.readlines()
    f.close()
    ias = []
    for line in lines:
        if line.startswith('#'): continue
        ias.append(line.split(',')[0])
    goods_composite_ia=0
    misses = 0
    ax3 = subplot(222)
    for ia in ias:
        idx = where(codex[:,1]==ia)
        if not codex[idx][:,0]: misses+=1; continue
        index = '%05d' %int(codex[idx][:,0][0])
        sfh = glob.glob('../SN_SFH/*/%s.dat'%index)
        if not sfh: continue
        data = loadtxt(sfh[0])[1:]
        ii = where((data[:,0]>lbl)&(data[:,0]<=lbu))
        data = data[ii]
        data[:,1] = boxcar(data[:,1],(bxsmooth,))
        goods_composite_ia += data[:,1]/simps(data[:,1])
        ax3.plot(data[:,0], data[:,1]/simps(data[:,1])/5., '-', color='0.65',alpha=0.3)
    norm1=len(ias)-misses
    ax3.plot(data[:,0],goods_composite_ia/norm1,'b-',label='%d GOODS SNe Ia' %norm1)
    ax3.plot(data[:,0],composite/len(sfhs),'k',label='All SNe')
    #test1 = goods_composite_ia/composite*(len(sfhs)/norm1)
    #ax3.plot(data[:,0], test1*0.005,'r-', label='GOODS Quotient')
    #ax3.axvspan(2.2,2.8,color='red',alpha=0.1)
    ax3.axvspan(1,4,color='red',alpha=0.1)
    

    f = open('candels_guesses.csv')
    lines = f.readlines()
    f.close()
    ias = []
    for line in lines:
        if line.startswith('#'): continue
        if (float(line.split(',')[-1])>0.5):
            ias.append(line.split(',')[0])
    candels_composite_ia=0
    misses = 0
    ax4 = subplot(223)
    for ia in ias:
        idx = where(codex[:,1]==ia)
        if not codex[idx][:,0]: misses+=1; continue
        index = '%05d' %int(codex[idx][:,0][0])
        sfh = glob.glob('../SN_SFH/*/%s.dat'%index)
        if not sfh: continue
        data = loadtxt(sfh[0])[1:]
        ii = where((data[:,0]>lbl)&(data[:,0]<=lbu))
        data = data[ii]
        data[:,1] = boxcar(data[:,1],(bxsmooth,))
        candels_composite_ia += data[:,1]/simps(data[:,1])
        ax4.plot(data[:,0], data[:,1]/simps(data[:,1])/5., '-', color='0.65',alpha=0.3)
    norm2=len(ias)-misses
    ax4.plot(data[:,0],candels_composite_ia/norm2,'b-',label='%d CANDELS SNe Ia' %norm2)
    ax4.plot(data[:,0],composite/len(sfhs),'k',label='All SNe')
    #test2 = candels_composite_ia/composite*(len(sfhs)/norm2)
    #ax4.plot(data[:,0], test2*0.005,'r-', label='CANDELS Quotient')
    ax4.axvspan(1,4,color='red',alpha=0.1)


    ax5 = subplot(224)
    ax5.plot(data[:,0],(candels_composite_ia+goods_composite_ia)/(norm1+norm2),
             'b-',label='%d SNe Ia' %(norm1+norm2))
    ax5.plot(data[:,0], composite/len(sfhs), 'k',label='All SNe')
    ax5.plot(data[:,0],(candels_composite_ia+goods_composite_ia)/composite*(len(sfhs)/(norm1+norm2))*0.005,
             'r-',label='SNe Ia Quotient')


    adjust_spines(ax1, ['bottom'])#['left'])
    adjust_spines(ax3, ['bottom'])
    adjust_spines(ax4, ['bottom'])#['left','bottom'])
    adjust_spines(ax5, ['bottom'])

    ylims = (0,0.01)
    xlims = (lbl,lbu)
    #ax1.set_ylim(ylims)
    #ax3.set_ylim(ylims)
    #ax4.set_ylim(ylims)
    #ax5.set_ylim(ylims)

    ax1.set_xlim(xlims)
    ax3.set_xlim(xlims)
    ax4.set_xlim(xlims)
    ax5.set_xlim(xlims)

    ax1.legend(loc=1,frameon=False)
    ax3.legend(loc=1,frameon=False)#, fontsize=10)
    ax4.legend(loc=1,frameon=False)#, fontsize=10)
    ax5.legend(loc=1,frameon=False)#, fontsize=10)
    

    plt.subplots_adjust(wspace=0.01,hspace=0.01)
    #show()
    savefig('plot1.png')#,transparent=True)
