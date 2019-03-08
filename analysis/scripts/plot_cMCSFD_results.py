#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from chainconsumer import ChainConsumer as chc
'''
For plotting the MCMC SFH results
L. Strolger
9/2018
'''



if __name__=='__main__':

    files = glob.glob('mc_sfd_*.pkl')
    tot = 0
    for file in files:
        samples = pickle.load(open(file,'rb'))
        ## samples = samples.reshape((-1,5))
        samples = samples[:,500:,:].reshape((-1,5))
        print('adding %d samples from %s... '%(len(samples), file))
        try:
            temp = concatenate((temp, samples), axis=0)
        except:
            temp = samples
        
    print('%d total samples' %len(temp))
    samples = temp
    

    ## gives back the 68% percentile range of values on each dimension
    ff_mcmc, m_mcmc, w_mcmc, k_mcmc, lnf_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                                    zip(*np.percentile(samples, [16, 50, 84],
                                                                       axis=0)))     
    print(r'ln scale: $\varepsilon=%2.3f\pm%2.3f$'%(ff_mcmc[0], ff_mcmc[1]))
    print(r'parameters: $\xi=%2.2f\pm%2.2f$; $\omega=%2.2f\pm%2.2f$; $\alpha=%2.2f\pm%2.2f$' %(m_mcmc[0], m_mcmc[1]
                                                                                               ,w_mcmc[0], w_mcmc[1]
                                                                                               ,k_mcmc[0], k_mcmc[1]))

    md0 = u.binmode(samples[:,1],bins=200)[0]
    md1 = u.binmode(samples[:,2],bins=200)[0]
    md2 = u.binmode(samples[:,3],bins=200)[0]
    md3 = u.binmode(samples[:,4],bins=200)[0]
    print(md0, md1, md2, md3)

    print()
    print('ff=%2.3f, %2.3f, %2.3f' %(ff_mcmc[0], ff_mcmc[1], ff_mcmc[2]))
    print('m=%2.2f, %2.2f, %2.2f' %(m_mcmc[0], m_mcmc[1], m_mcmc[2]))
    print('w=%2.2f, %2.2f, %2.2f' %(w_mcmc[0], w_mcmc[1], w_mcmc[2]))
    print('k=%2.2f, %2.2f, %2.2f' %(k_mcmc[0], k_mcmc[1], k_mcmc[2]))
    print('lnf=%2.2f, %2.2f, %2.2f' %(lnf_mcmc[0], lnf_mcmc[1], lnf_mcmc[2]))

    parameters = [r'$\log\varepsilon$',r'$\xi$',r'$\omega$',r'$\alpha$', r'$\log f$']
    c = chc()
    c.add_chain(samples, parameters=parameters, walkers=1000, name='Unimodal Model')

    ## table1 = c.analysis.get_latex_table(caption="Results for the tested model", label="tab:example")
    ## print(table1)

    ## latex_table = c.analysis.get_correlation_table()
    ## print (latex_table)

    ## latex_table = c.analysis.get_covariance_table()
    ## print(latex_table)

    ## fig = c.plotter.plot_walks(truth={r'$\varepsilon$':ff_mcmc[0], 
    ##                                   r'$\xi$': md0, r'$\omega$': md1, r'$\alpha$': md2,
    ##                                   r'$\log f$':lnf_mcmc[0]},
    ##                            convolve=100)
    ## ## fig = c.plotter.plot_walks(parameters=parameters[1:4], convolve=1000)
    ## savefig('temp.png')
    ## sys.exit()

    ## fig = c.plotter.plot_distributions(truth={r'$\xi$': md0, r'$\omega$': md1, '$\alpha$': k_mcmc[0]})

    ## grc = c.diagnostic.gelman_rubin()
    ## gc = c.diagnostic.geweke()
    ## print(grc,gc)

    c.configure(label_font_size=18, contour_labels='sigma')
    fig = c.plotter.plot(figsize="column",
                         truth={r'$\log\varepsilon$':ff_mcmc[0], 
                                r'$\xi$': m_mcmc[0], r'$\omega$': w_mcmc[0], r'$\alpha$': k_mcmc[0],
                                r'$\log f$':md3},
                         extents = [[-3.1, -2.5], [-2000.,0.0], [10., 90], [0., 500],[-4,0]]
                         ## extents = [[-0.2, 0.4], [-1900.,200.0], [-0.1, 90], [-200., 500],[-4.2,0.1]]
                         )
    fig.set_size_inches(4.5 + fig.get_size_inches())
    savefig('figure_sfd_corners.png')
    
