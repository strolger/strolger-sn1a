#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import imf
from scipy.integrate import quad


def imfs(m1,m2):
    num1 = quad(imf.salpeter,m1,m2)[0]
    den1 = quad(imf.salpeter1,0.1,125)[0]
    ## print ('Salpeter k=%1.4f' %(num1/den1))
    

    p0=[0.5,1.]
    num = quad(imf.kroupa,m1,m2,args=tuple(p0))[0]
    den = quad(imf.kroupa1,0.1,125,args=tuple(p0))[0]
    ## print ('Kroupa k=%1.4f' %(num/den))
    return(num1/den1, num/den)


data = []
for m1 in arange(2.0,3.5,0.025):
    for m2 in arange(6,11.0,0.025):
        data.append(imfs(m1,m2))

data=array(data)
med=quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
print(med)
#print(average(data[:,0]), average(data[:,1]))
#print(std(data[:,0]), std(data[:,1]))
#print('minus %2.1f%%, plus %2.1f%%' %(std(data[:,0])/med*100., std(data[:,1])/med*100.))
bins =100
N1, bins=histogram(data[:,0],bins=bins)
N2, bins=histogram(data[:,1],bins=bins)

N1 = N1/sum(N1)
N2 = N2/sum(N2)

ax=subplot(111)
ax.plot(bins[:-1],0.5*(N1+N2), 'ko')

yy,xx = zip(*sorted(zip((N1+N2)/2.,bins[:-1])))
yy = array(yy)[::-1]
xx = array(xx)[::-1]
intv=0.68

idx = where(cumsum(yy)<=intv)
print(min(xx[idx]), max(xx[idx]))
ax.axvline(med)
ax.axvline(min(xx[idx]),color='red')
ax.axvline(max(xx[idx]),color='red')
print('minus %2.1f%%, plus %2.1f%%' %((med-min(xx[idx]))/med*100., (max(xx[idx])-med)/med*100.))
idx = where(yy.max())
ax.axvline(xx[idx],linestyle='--', color='blue')
print(xx[idx])
ax.set_xlabel('k-value')
ax.set_ylabel('Normalized distribution')
ax.set_xlim(0,0.05)
savefig('junk.png')

pdb.set_trace()
