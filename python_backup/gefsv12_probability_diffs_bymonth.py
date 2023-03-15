"""
python gefsv12_probability_diffs_bymonth.py cyyyymmddhh clead cthresh
"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import _pickle as cPickle
import scipy.stats as stats
from dateutils import dateshift, daterange

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2] # 018 not 18
cthresh = sys.argv[3]  #0.254, 1.0,  5.0, 10.0, 25.0
jnd = 1350
ind = 1110

cyear = cyyyymmddhh[0:4]
rthresh = float(cthresh)
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ndaysomo = [31,28,31,30,31,30, 31,31,30,31,30,31]
ccmonth = cmonths[imonth]
cdomain = 'conus'
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/all_thinned/'
cyyyymmddhh_begin = cyyyymmddhh
ndays = ndaysomo[imonth]-1
cyyyymmddhh_end = dateshift(cyyyymmddhh_begin, ndays*24)
cyyyymmddhh_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)

# ---- read the stored netCDF thinned probabilities

infile = master_directory_thinned_output +\
    ccmonth + cyear+'_lead'+clead+\
            '_probabilities_all_thinned.nc'
print ('reading from ', infile)
nc = Dataset(infile,'r')
yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
print ('yyyymmddhh = ', yyyymmddhh)
thresholds = nc.variables['thresholdv'][:]
if rthresh == 0.254:
    itt = 0
else:
    itt = int(np.where(thresholds == rthresh)[0])
print ('itt = ', itt)
lons_ndfd = nc.variables['lons'][:,:]
lats_ndfd = nc.variables['lats'][:,:]
for iday, cday in enumerate(cyyyymmddhh_list):
    idd = np.where(yyyymmddhh == int(cday))[0]
    probability_qmapped_weighted = np.squeeze(\
        nc.variables['probability_qmapped_weighted'][idd,itt,:,:])
    probability_qmapped_weighted_dressed = np.squeeze(\
        nc.variables['probability_qmapped_weighted_dressed'][idd,itt,:,:])
    diff = probability_qmapped_weighted_dressed - probability_qmapped_weighted
    print ('cday, max, min diff = ', cday, np.max(diff), np.min(diff))

nc.close()
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)

