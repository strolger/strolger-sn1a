#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
import corner
'''
For plotting the MCMC SFH results
L. Strolger
5/2018
'''



if __name__=='__main__':

    ## something to work the output of hte sampler.
    ''' use corner, and play.
    
    https://corner.readthedocs.io/en/latest/pages/custom.html
    '''

    files = glob.glob('mc_sfh_201806*.pkl')
    temp = zeros((1,3),)
    tot = 0
    for file in files:
        samples = pickle.load(open(file,'rb'))
        ## ax1 = subplot(311)
        ## ax2 = subplot(312)
        ## ax3 = subplot(313)
        ## for i in range(len(samples[:,:,0])):
        ##     ax1.plot(range(len(samples[:,:,0][0,:])), samples[:,:,0][i,:],'k-')
        ##     ax2.plot(range(len(samples[:,:,0][0,:])), samples[:,:,1][i,:],'k-')
        ##     ax3.plot(range(len(samples[:,:,0][0,:])), samples[:,:,2][i,:],'k-')
        ## show()
        ## clf()
        samples = samples[:,150:,:].reshape((-1,3))
        print('adding %d samples from %s... '%(len(samples), file))
        temp = concatenate((temp, samples), axis=0)
        
    print('%d total samples' %len(temp))

    samples = temp

    ## gives back the 68% percentile range of values on each dimension
    m_mcmc, w_mcmc, k_mcmc = map(lambda v: (v[1], v[1]-v[0], v[2]-v[1]),
                                 zip(*np.percentile(samples, [16, 50, 84],
                                                    axis=0)))     
    ## print(r'parameters: $\xi=%2.2f\pm%2.2f$; $\omega=%2.2f\pm%2.2f$; $\alpha=%2.2f\pm%2.2f$' %(m_mcmc[0], m_mcmc[1]
    ##                                                                                            ,w_mcmc[0], w_mcmc[1]
    ##                                                                                            ,k_mcmc[0], k_mcmc[1]))

    md0 = u.binmode(samples[:,0],bins=20)[0]
    md1 = u.binmode(samples[:,1],bins=20)[0]
    md2 = u.binmode(samples[:,2],bins=20)[0]

    print(md0, md1, md2)

    print('m=%2.2f, %2.2f, %2.2f' %(m_mcmc[0], m_mcmc[1], m_mcmc[2]))
    print('w=%2.2f, %2.2f, %2.2f' %(w_mcmc[0], w_mcmc[1], w_mcmc[2]))
    print('k=%2.2f, %2.2f, %2.2f' %(k_mcmc[0], k_mcmc[1], k_mcmc[2]))

    import corner
    fig = corner.corner(samples,labels=[r'$\xi$',r'$\omega$',r'$\alpha$'],
                        truths=[md0, md1, md2])
    

    
                    

    ## value1 = array([md0,md1,md2])
    ## value2 = samples.mean(axis=0)
    ## ndim = 3
    ## # Extract the axes
    ## axes = np.array(fig.axes).reshape((ndim, ndim))

    ## # Loop over the diagonal
    ## for i in range(ndim):
    ##     ax = axes[i, i]
    ##     ax.axvline(value1[i], color="g")
    ##     ax.axvline(value2[i], color="r")

    ##     # Loop over the histograms
    ##     for yi in range(ndim):
    ##         for xi in range(yi):
    ##             ax = axes[yi, xi]
    ##             ax.axvline(value1[xi], color="g")
    ##             ax.axvline(value2[xi], color="r")
    ##             ax.axhline(value1[yi], color="g")
    ##             ax.axhline(value2[yi], color="r")
    ##             ax.plot(value1[xi], value1[yi], "sg")
    ##             ax.plot(value2[xi], value2[yi], "sr")    

    show()
