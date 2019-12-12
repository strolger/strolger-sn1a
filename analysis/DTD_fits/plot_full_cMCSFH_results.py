#!/usr/bin/env python
import os,sys,pdb,scipy,glob,pickle
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz
from chainconsumer import ChainConsumer as chc
'''
For plotting the temp MCMC SFH results
L. Strolger
12/2019
'''



if __name__=='__main__':

    ## files = glob.glob('mc_sfh_201810152324.pkl')
    #files = glob.glob('mc_sfh_2019*.pkl')
    #files = ['mc_sfh_201912061430.pkl']
    #files = ['mc_sfh_201912091319.pkl']
    files = glob.glob('mc_sfh_201912*.pkl')
    tot = 0
    for file in files:
        samples = pickle.load(open(file,'rb'))
        samples = samples[:,50:,:].reshape((-1,4))
        ## samples = samples.reshape((-1,3))
        print('adding %d samples from %s... '%(len(samples), file))
        try:
            temp = concatenate((temp, samples), axis=0)
        except:
            temp = samples
        
    print('%d total samples' %len(temp))
    samples = temp

    ## gives back the 68% percentile range of values on each dimension
    f_mcmc, m_mcmc, w_mcmc, k_mcmc = map(lambda v: (v[1], v[1]-v[0], v[2]-v[1]),
                                         zip(*np.percentile(samples, [16, 50, 84],
                                                            axis=0)))     
    print(r'parameters: $\ln\varepsilon=%2.2f\pm%2.2f$; $\xi=%2.2f\pm%2.2f$; $\omega=%2.2f\pm%2.2f$; $\alpha=%2.2f\pm%2.2f$' %(f_mcmc[0], f_mcmc[1]
                                                                                                                ,m_mcmc[0], m_mcmc[1]
                                                                                                                ,w_mcmc[0], w_mcmc[1]
                                                                                                                ,k_mcmc[0], k_mcmc[1]))
    
    md0 = u.binmode(samples[:,0],bins=20)[0]
    md1 = u.binmode(samples[:,1],bins=20)[0]
    md2 = u.binmode(samples[:,2],bins=20)[0]
    md3 = u.binmode(samples[:,2],bins=20)[0]

    print(md0, md1, md2, md3)

    print('f=%2.2f, %2.2f, %2.2f' %(f_mcmc[0], f_mcmc[1], f_mcmc[2]))
    print('m=%2.2f, %2.2f, %2.2f' %(m_mcmc[0], m_mcmc[1], m_mcmc[2]))
    print('w=%2.2f, %2.2f, %2.2f' %(w_mcmc[0], w_mcmc[1], w_mcmc[2]))
    print('k=%2.2f, %2.2f, %2.2f' %(k_mcmc[0], k_mcmc[1], k_mcmc[2]))

    parameters = [r'$\ln\varepsilon$',r'$\xi$',r'$\omega$',r'$\alpha$']
    c = chc()
    c.add_chain(samples, parameters=parameters, walkers=100, name='SFH MCMC')

    latex_table = c.analysis.get_correlation_table()
    print (latex_table)
    #latex_table = c.analysis.get_covariance_table()
    #print(latex_table)
    #pdb.set_trace()
    ## fig = c.plotter.plot_walks(parameters=parameters[:1], convolve=100)
    ## fig = c.plotter.plot_walks(truth={r'$\xi$': m_mcmc[0], r'$\omega$': w_mcmc[0], '$\alpha$': w_mcmc[0]},
    ##                           convolve=1000)
    ## savefig('temp.png')
    ## sys.exit()
    ## fig = c.plotter.plot_distributions(truth={r'$\xi$': md0, r'$\omega$': md1, '$\alpha$': k_mcmc[0]})

    ## grc = c.diagnostic.gelman_rubin()
    ## gc = c.diagnostic.geweke()
    ## print(grc,gc)


    
    marker_data = array([[-2.78, -1518, 51, 50],
                         [-2.78, -1518, 51, 50]])
        
    c.add_chain(marker_data, posterior=np.array([0,1]), plot_contour=False, 
                plot_point=True, marker_size=170, marker_style="o", color="r",
                name="CSFH Max. Likelihood (Optimized)")

    c.configure(label_font_size=22, tick_font_size=14, contour_labels='sigma')
    fig = c.plotter.plot(figsize="column",
                         truth={r'$\ln\varepsilon$':md0, r'$\xi$': md1, r'$\omega$': md2, '$\alpha$': k_mcmc[0]},
                         extents = [[-3.5,-2.5],[-1900.,-100.0], [-10., 90], [-10., 500]]
                         )
    fig.set_size_inches(4.5 + fig.get_size_inches())
    savefig('figure_full_sfh_corners.pdf')
    
