#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from scipy.stats import gaussian_kde
from random import sample


## lets define some plots

rcParams['figure.figsize']=20, 10
rcParams['font.size']=24.0

fig = plt.figure()
gs = GridSpec(4,4)

ax1 = fig.add_subplot(gs[1:4,0:3])
ax2 = fig.add_subplot(gs[0,0:3], sharex=ax1)
ax3 = fig.add_subplot(gs[1:4,3], sharey=ax1)



#loadin up data
file1 = '../ALLSFH_new_z/Cami_GOODS-N_zbest.dat'
file2 = '../ALLSFH_new_z/Cami_GOODS-S_zbest.dat'

data1 = loadtxt(file1)
data1[:,0] = data1[:,0]+100000
data2 = loadtxt(file2)
data2[:,0] = data2[:,0]+200000
data = concatenate((data1,data2), axis=0)

## getting rid of bad data
idx = where(data[:,1]==-99)
data[idx[0]]=float('nan')
idx = where(data[:,1]==8.0)
data[idx[0]]=float('nan')
idx = where(data[:,3]==-99)
data[idx[0]]=float('nan')
idx = where(data[:,3]==-3.0)
data[idx[0]]=float('nan')

## and to get rid of stars
"""
This is tricky as the info for CLASS_STAR
is in a separate set of files
"""

j1 = loadtxt('../ALLSFH_new_z/CANDELS_GDSN_znew_avgal_radec.dat')
j2 = loadtxt('../ALLSFH_new_z/CANDELS_GDSS_znew_avgal_radec.dat')
j1 = np.delete(j1,[38,39],1) ## remove the extra phot columns 

## j3 = zeros((len(j2), len(j1[0])),)
## j3[:,:len(j2[0])]=j2
## j3[:,-3:]=j2[:,-3:]
           

j1[:,0] = j1[:,0]+100000
j2[:,0] = j2[:,0]+200000
junk = concatenate((j1,j2), axis=0)
idx = where(junk[:,-3]>=0.8)
data[idx[0]]=float('nan')



X = data[:,1][~isnan(data[:,1])]
Y = data[:,3][~isnan(data[:,3])]
all_masses = X
all_sfrs = Y
print(len(X), len(Y))
#sys.exit()

# fit an array of size [Ndim, Nsamples]
tmp = np.vstack([X, Y])
kde = gaussian_kde(tmp)


# evaluate on a regular grid
xgrid = np.linspace(8.0, 12.,60)
ygrid = np.linspace(-3,4,60) 
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))


ax1.imshow(Z.reshape(Xgrid.shape),
           origin='lower', aspect='auto',
           extent=[8,12, -3, 4],
           ##cmap=plt.cm.gist_earth_r)#, alpha=0.6)
           cmap='bone_r')
#cb = plt.colorbar()
#cb.set_label("density")

cs1 = ax1.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape), [0.01, 0.05, 0.32], colors='black')
fmt={}
strs = ['3$\sigma$','2$\sigma$','1$\sigma$']
for l, s in zip(cs1.levels, strs):
    fmt[l]=s
plt.clabel(cs1, cs1.levels, inline=1, fmt=fmt)

## ii = range(len(X))
## ix = sample(ii,1000)
## plt.plot(X[ix],Y[ix], 'ko')

ax2.hist(all_masses, bins=100,  color='0.5', normed=True, label='All Catalog Galaxies')
ax3.hist(all_sfrs, bins=100,  color='0.5', normed=True, orientation='horizontal')

idx = loadtxt('../DTD_fits/host_idxs.txt')
mask = in1d(data[:,0], idx)
data_ias = data[mask]


X = data_ias[:,1]
Y = data_ias[:,3]
host_masses = X
host_sfrs = Y

# fit an array of size [Ndim, Nsamples]
tmp = np.vstack([X, Y])
kde = gaussian_kde(tmp)


# evaluate on a regular grid
xgrid = np.linspace(8.0, 12.,40)
ygrid = np.linspace(-3,4,40) 
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))

## # Plot the result as an image
## plt.imshow(Z.reshape(Xgrid.shape),
##            origin='lower', aspect='auto',
##            extent=[8,12, -3, 4],
##            ## cmap=plt.cm.gist_earth_r, alpha=0.4)
##            cmap='Blues', alpha=0.5)
## ## cb = plt.colorbar()
## ## cb.set_label("density")
## ## pdb.set_trace()

#cs2=ax1.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape),[0.01, 0.05, 0.2, 0.32],  colors='C3', linestyle='dashed')
## cs=plt.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape),[0.01, 0.05, 0.32],  colors='blue')
## fmt={}
## strs = ['3$\sigma$','2$\sigma$','1$\sigma$']
## for l, s in zip(cs.levels, strs):
##     fmt[l]=s
## plt.clabel(cs, cs.levels, inline=1, fmt=fmt)

ax1.plot(X,Y, 'o', color='C3', ms=10)
ax2.hist(host_masses, bins=10,  align='left', color ='C3', lw=3, histtype = 'step', normed=True, label='SN Ia Hosts')
ax3.hist(host_sfrs, bins=10,  color='C3', lw=3, histtype= 'step', normed=True, orientation='horizontal')


ax2.legend(loc=1, ncol=2, frameon=False, fontsize=20)

ax1.set_xlim(8,12)
ax2.set_ylim(0,0.8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
ax2.set_yticks([])
ax3.set_xticks([])

ax1.set_xlabel('Log (M/M$_{\odot}$)')
ax1.set_ylabel('Log (SFR)')
savefig('figure_contours.pdf')
