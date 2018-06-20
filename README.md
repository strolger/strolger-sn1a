# strolger-sn1a

This analysis looks to find a maximum likelihood delay time distribution
for type Ia supernovae from the convolution of DTD models with: (a) star
formation rate density histories to match volumetric SN Ia rates, and
(b) star formation histories of individual galaxies match expected
number of events to the observed number over the duration of a survey.

There are coded methods to find the former through both optimized and
MCMC methods, and the latter through MCMC.


And now, a brief tour...

README.md-- (this file) an introduction to the analysis

notes.txt-- an explaination of what this analysis is designed to test
and some caveats.

paper/-- will be the formal writeup. There are some scratches there.

analysis/-- is where the meat of this work is.

analysis/older_SN_SFHs.tgz-- Some preliminary SFHs C. Pacifici provided
to draft code from

analysis/ALLSFH_new_z-- contains the SFHs and catalog information
C. Pacifici provided.

analysis/simple_plots-- some preliminary plots of the SN SFHs in
comparison to the population of galaxies.

analysis/scripts-- where the REAL action is!
There's a lot of scripts there, but for now I'll just focus on the main
MCMC... MCMC_SFH.py. This does is the test of the DTD on the SFHs. It's
somewhat self explanitory. Two outputs are temporary.png, which shows a
very preliminary corner plot. The other is a pickle file of the samples,
prefixed by "mc_sfh_" and encoded with the date-time of the run. 

plot_cMCSFH_results.py gives a slightly more detailed corner plot,
figure_sfh_corners.png, and other analysis. As a note, it picks up all
mcmc pickle files.


contact me: strolger@stsci.edu
