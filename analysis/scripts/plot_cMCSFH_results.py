#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from chainconsumer import ChainConsumer as chc
'''
For plotting the MCMC SFH results
L. Strolger
5/2018
'''



if __name__=='__main__':

    files = glob.glob('mc_sfh_201806*.pkl')
    ## temp = zeros((1,3),)
    tot = 0
    for file in files:
        samples = pickle.load(open(file,'rb'))
        samples = samples[:,100:,:].reshape((-1,3))
        print('adding %d samples from %s... '%(len(samples), file))
        try:
            temp = concatenate((temp, samples), axis=0)
        except:
            temp = samples
        
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

    parameters = [r'$\xi$',r'$\omega$',r'$\alpha$']
    c = chc()
    c.add_chain(samples, parameters=parameters, walkers=60)
    latex_table = c.analysis.get_correlation_table()
    print (latex_table)
    latex_table = c.analysis.get_covariance_table()
    print(latex_table)
    ## fig = c.plotter.plot_walks(truth={r'$\xi$': md0, r'$\omega$': md1, '$\alpha$': md2},
    ##                            convolve=100)
    ## show()
    ## sys.exit()
    ## fig = c.plotter.plot_distributions(truth={r'$\xi$': md0, r'$\omega$': md1, '$\alpha$': k_mcmc[0]})

    ## grc = c.diagnostic.gelman_rubin()
    ## gc = c.diagnostic.geweke()
    ## print(grc,gc)

    fig = c.plotter.plot(figsize="column", truth={r'$\xi$': md0, r'$\omega$': md1, r'$\alpha$': md2})
    fig.set_size_inches(4.5 + fig.get_size_inches())
    savefig('figure_sfh_corners.png')
    
