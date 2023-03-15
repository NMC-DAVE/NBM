"""
decayavg_filter_conus_B_2018.py
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
                
# ---- read covariance matrices from cPickle file

#infile = 'covstats_bias_random_ecmwf2019_localized_lead='+clead+'.cPick'
infile = 'covstats_bias_random_ecmwf2018_localized_lead='+clead+'.cPick'
inf = open(infile,'rb')
cov_bias_localized = cPickle.load(inf) 
cov_random_localized = cPickle.load(inf) 
var_bias = cPickle.load(inf) 
var_random = cPickle.load(inf) 
lats = cPickle.load(inf) 
lons = cPickle.load(inf) 
nlats, nlons = np.shape(lats)
inf.close()  

# --- make the sum of the random and the bias into a 2D-array.
#     to the diagonal add the obs variance.   Then invert

npts = nlats*nlons
B_and_R = np.zeros((npts, npts), dtype=np.float64)
ktr1 = 0
for i1 in range(nlons):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('processing i1 = ',i1,' of ',nlons,'. Current time = ',current_time)
    for j1 in range(nlats):
        ktr2 = 0 
        for i2 in range(nlons): 
            for j2 in range(nlats):
                B_and_R[ktr1,ktr2] = cov_bias_localized[j1,i1,j2,i2] + \
                    1.5*cov_random_localized[j1,i1,j2,i2] 
                ktr2 = ktr2 + 1   
        ktr1 = ktr1 + 1
        
#for ktr1 in range(npts):
#    B_and_R[ktr1,ktr1] = B_and_R[ktr1,ktr1] + 1.0

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('Inverting matrix. Current time = ',current_time)
B_and_R_inverse = np.linalg.inv(B_and_R) 
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('Finished inverting matrix. Current time = ',current_time)

# ---- reform matrix inverse into 4D-array.

matrix_inverse_4D = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float64)
ktr1 = 0
for i1 in range(nlons):
    for j1 in range(nlats):
        ktr2 = 0 
        for i2 in range(nlons): 
            for j2 in range(nlats):
                matrix_inverse_4D[j1,i1,j2,i2] = B_and_R_inverse[ktr1,ktr2] 
                ktr2 = ktr2+1   
        ktr1 = ktr1 + 1
        
# ---- save to cPickle file.

#outfile = 'B_and_R_matrix_inverse_lead='+clead+'.cPick'
outfile = 'B_and_R_matrix_inverse_2018_lead='+clead+'.cPick'
print (outfile)
ouf = open(outfile, 'wb')
cPickle.dump(matrix_inverse_4D, ouf)
ouf.close()

    