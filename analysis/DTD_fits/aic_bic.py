#!/usr/bin/env python
'''
Simple AIC/BIC
L. Strolger
2019
'''
import os,sys,pdb,scipy,glob,pickle,datetime
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from strolger_util import imf
from strolger_util import cosmotools as ct
from copy import copy, deepcopy
from scipy.integrate import simps,quad
import emcee
import warnings
warnings.simplefilter("ignore",RuntimeWarning)



def aic_bic(x,y,ym,nvar):
    ssc = sum((y - ym)**2.)
    aic = 2.*nvar-2.*log(ssc)
    bic = len(y) * log(ssc/len(y)) + nvar*log(ssc)
    return(aic,bic)

if __name__=='__main__':

    calculate = False
    
    rates = loadtxt('SNeIa_rates.txt')
    rates[:,1:] = rates[:,1:]#*1.0e-4 ## put on the right scale
    rates = rates[:,:4]
    brates = u.gimme_rebinned_data(rates,verbose=False,splits=arange(0,1.167,0.167).tolist())
    data = deepcopy(rates) #this seems to work better than then actual values

    p0 = [-1518, 51,  50]
    p1 = [-1258, 59, 248]
    p2 = [-1087, 73, 242]
    frac = 0.062

    age = 13.6
    scale_k = quad(imf.salpeter,3,8)[0]/quad(imf.salpeter1,0.1,125)[0]
    scale = scale_k * 0.7**2.*1e4## factors of h...
    dt = 0.05
    tt = arange(dt,age,dt)
    lbt = age - tt
    zz = [ct.cosmoz(x, ho=70) for x in lbt]

    par_model = [0.0134, 2.55, 3.3, 6.1]

    sfh = rz.csfh_time(tt, *par_model)
    pwrl = (-1.0,1.0)
    dud = rz.powerdtd(tt, *pwrl)

    dtd = rz.dtdfunc(tt, *p0)
    tmp = convolve(sfh, dtd, 'full')*dt*frac*scale
    rate_fn=tmp[:len(dtd)]
    junk, m_rates = u.recast(rates[:,0], 0., zz, rate_fn) ## optimized fit
    (aic, bic)= aic_bic(rates[:,0], rates[:,1], m_rates, 4)
    print('Optimized fit: aic=%2.1f, bic=%2.1f'%(aic,bic))


    jdud = convolve(sfh, dud, 'full')*dt*scale*0.065
    jdud = jdud[:len(dtd)]
    junk, m_rates = u.recast(rates[:,0], 0., zz, jdud) ## power law
    (aic, bic)= aic_bic(rates[:,0], rates[:,1], m_rates, 2)
    print('power-law fit: aic=%2.1f, bic=%2.1f'%(aic,bic))



    dtd = rz.dtdfunc(tt, *p1)
    tmp = convolve(sfh, dtd, 'full')*dt*frac*scale
    rate_fn=tmp[:len(dtd)]
    junk, m_rates = u.recast(rates[:,0], 0., zz, rate_fn) ## median CSFH
    (aic, bic)= aic_bic(rates[:,0], rates[:,1], m_rates, 4)
    print('Median CSFH fit: aic=%2.1f, bic=%2.1f'%(aic,bic))

    dtd = rz.dtdfunc(tt, *p2)
    tmp = convolve(sfh, dtd, 'full')*dt*frac*scale
    rate_fn=tmp[:len(dtd)]
    junk, m_rates = u.recast(rates[:,0], 0., zz, rate_fn) ## median SFH
    (aic, bic)= aic_bic(rates[:,0], rates[:,1], m_rates, 4)
    print('Median SFH fit: aic=%2.1f, bic=%2.1f'%(aic,bic))
