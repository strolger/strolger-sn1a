#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from scipy.stats import gaussian_kde
from random import sample


rcParams['figure.figsize']=20, 10
rcParams['font.size']=24.0

## all_mass = loadtxt('all_masses.txt')
## host_mass = loadtxt('host_masses.txt')
## all_sfr = loadtxt('all_sfrs.txt')
## host_sfr = loadtxt('host_sfrs.txt')


file1 = '../ALLSFH_new_z/Cami_GOODS-N_zbest.dat'
file2 = '../ALLSFH_new_z/Cami_GOODS-S_zbest.dat'

data1 = loadtxt(file1)
data1[:,0] = data1[:,0]+100000
data2 = loadtxt(file2)
data2[:,0] = data2[:,0]+200000
data = concatenate((data1,data2), axis=0)



idx = where(data[:,1]==-99)
data[idx[0]]=float('nan')
idx = where(data[:,1]==8.0)
data[idx[0]]=float('nan')
idx = where(data[:,3]==-99)
data[idx[0]]=float('nan')
idx = where(data[:,3]==-3.0)
data[idx[0]]=float('nan')

X = data[:,1][~isnan(data[:,1])]
Y = data[:,3][~isnan(data[:,3])]

# fit an array of size [Ndim, Nsamples]
tmp = np.vstack([X, Y])
kde = gaussian_kde(tmp)


# evaluate on a regular grid
xgrid = np.linspace(8.0, 12.,40)
ygrid = np.linspace(-3,4,40) 
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))


# Plot the result as an image
plt.imshow(Z.reshape(Xgrid.shape),
           origin='lower', aspect='auto',
           extent=[8,12, -3, 4],
           ##cmap=plt.cm.gist_earth_r)#, alpha=0.6)
           cmap='bone_r')
#cb = plt.colorbar()
#cb.set_label("density")
## pdb.set_trace()

cs = plt.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape), [0.01, 0.05, 0.32], colors='black')
## fmt={}
## strs = ['3$\sigma$','2$\sigma$','1$\sigma$']
## for l, s in zip(cs.levels, strs):
##     fmt[l]=s
## plt.clabel(cs, cs.levels, inline=1, fmt=fmt)

## ii = range(len(X))
## ix = sample(ii,1000)
## plt.plot(X[ix],Y[ix], 'ko')



idx = loadtxt('../DTD_fits/host_idxs.txt')
mask = in1d(data[:,0], idx)
data_ias = data[mask]


X = data_ias[:,1]
Y = data_ias[:,3]

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

cs=plt.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape),[0.01, 0.05, 0.2, 0.32],  colors='blue', linestyle='dashed')
## cs=plt.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape),[0.01, 0.05, 0.32],  colors='blue')
## fmt={}
## strs = ['3$\sigma$','2$\sigma$','1$\sigma$']
## for l, s in zip(cs.levels, strs):
##     fmt[l]=s
## plt.clabel(cs, cs.levels, inline=1, fmt=fmt)

#plt.plot(X,Y, 'bo')
savefig('figure_contor.png')
