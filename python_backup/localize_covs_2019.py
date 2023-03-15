"""

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.signal as signal
import scipy.stats as stats
from astropy.convolution import convolve


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

clead = sys.argv[1]
cefolda = sys.argv[2] # e-folding distance in grid units.
cefoldb = sys.argv[3]
cefoldg = sys.argv[4]
efolda = float(cefolda)
efoldb = float(cefoldb)
efoldg = float(cefoldg)
efold2a = float(efolda)*float(efolda)
efold2b = float(efoldb)*float(efoldb)
                
# ---- read from cPickle file

#infile = 'covstats_bias_random_ecmwf2019_lead='+clead+'.cPick'
infile = 'covstats_bias_random_ecmwf2018_lead='+clead+'.cPick'
inf = open(infile,'rb')
cov_bias = cPickle.load(inf)
cov_random = cPickle.load(inf)
var_bias = cPickle.load(inf)
var_random = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
nlats, nlons = np.shape(lats)
inf.close()         

# ---- localize

degtorad = 3.1415926/180.
#gaussian_b = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
#gaussian_r = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
exponential_b = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
exponential_r = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
analysis_error_cov = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
for ilat2 in range(nlats):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('processing ilat2 = ',ilat2,' of ',nlats,'. Current time = ',current_time)
    for ilon2 in range(nlons):
        for ilat1 in range(nlats):
            for ilon1 in range(nlons):
                idiff = np.abs(lons[ilat2,ilon2]-lons[ilat1,ilon1])
                jdiff = np.abs(lats[ilat2,ilon2]-lats[ilat1,ilon1])
                coslat = np.cos(lats[ilat2,ilon2]*degtorad) + \
                    np.cos(lats[ilat1,ilon1]*degtorad)
                #dist2 = (idiff*coslat)**2 + jdiff**2
                dist = np.sqrt((idiff*coslat)**2 + jdiff**2)
                #gaussian_b[ilat1,ilon1,ilat2,ilon2] = \
                #    np.exp(-dist2/efold2a)
                #gaussian_r[ilat1,ilon1,ilat2,ilon2] = \
                #    np.exp(-dist2/efold2b)
                exponential_b[ilat1,ilon1,ilat2,ilon2] = np.exp(-dist/efolda)
                exponential_r[ilat1,ilon1,ilat2,ilon2] = np.exp(-dist/efoldb)
                analysis_error_cov[ilat1,ilon1,ilat2,ilon2] = np.exp(-dist**2/efoldg**2)

# --- product of two

#cov_bias_localized = cov_bias*gaussian_b
#cov_random_localized = cov_random*gaussian_r 
cov_bias_localized = cov_bias*exponential_b
cov_random_localized = cov_random*exponential_r 

# --- write out

#outfile = 'covstats_bias_random_ecmwf2019_localized_lead='+clead+'.cPick'
outfile = 'covstats_bias_random_ecmwf2019_localized_lead='+clead+'.cPick'
ouf = open(outfile,'wb')
cPickle.dump(cov_bias_localized, ouf)
cPickle.dump(cov_random_localized, ouf)
cPickle.dump(analysis_error_cov, ouf)
cPickle.dump(var_bias, ouf)
cPickle.dump(var_random, ouf)
cPickle.dump(lats, ouf)
cPickle.dump(lons, ouf)
ouf.close()       
