# bayesianblocks

bayesian_lc is a user-oriented python 3 program designed to construct bayesian blocks [1] — optimal piecewise-constant approximation of time series. Originally written for operating with Fermi/LAT 4FGL-DR2 light curves [2], the script might be easily edited for various datasets. The principal part of the program consists of the optimal segmentation [3] and the function block_av which averages the flux in each bin. Plotting is also implemented.

[1] https://ui.adsabs.harvard.edu/abs/2013ApJ...764..167S

[2] https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

[3] https://docs.astropy.org/en/stable/_modules/astropy/stats/bayesian_blocks.html
