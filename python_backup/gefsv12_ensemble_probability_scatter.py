"""
python gefsv12_ensemble_probability_weighted.py cyyyymmddhh clead cthresh
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
ccmonth = cmonths[imonth]
cdomain = 'conus'
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'

plot_fullfield = True

# ---- read the stored netCDF full-field probabilities, quantile mapped

if plot_fullfield == True:
    infile = master_directory_fullfield_output +\
        ccmonth+cyear+'_use99_lead'+clead+\
        '_probabilities_weighted_fullfield.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
    idd = np.where(yyyymmddhh == int(cyyyymmddhh))[0]
    print ('yyyymmddhh = ', yyyymmddhh)
    thresholds = nc.variables['thresholdv'][:]
    if rthresh == 0.254:
        itt = 0
    else:
        itt = int(np.where(thresholds == rthresh)[0])
    print ('itt = ', itt)
    lons_ndfd = nc.variables['lons'][:,:]
    lats_ndfd = nc.variables['lats'][:,:]
    prob_qmapped_fullfield = \
        np.squeeze(nc.variables['probability_qmapped'][idd,itt,:,:])
    nc.close()
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    
# ---- read the stored netCDF full-field probabilities, quantile mapped

if plot_fullfield == True:
    infile = master_directory_fullfield_output +\
        ccmonth+cyear+'_use99_lead'+clead+\
        '_probabilities_weighted_dressed_fullfield.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    prob_qmapped_fullfield_dressed = \
        np.squeeze(nc.variables['probability_qmapped'][idd,itt,:,:])
    nc.close()

# ---- make scatterplot of probabilities

fig = plt.figure(figsize=(6.,6.5))
axloc = [0.15,0.15,0.8,0.76]
ax = fig.add_axes(axloc)
cleadb = str(int(clead)-6)

title = 'Scatterplot, QM weighted & dressed vs. QM weighted,\n'+cyyyymmddhh+', lead = '+clead+' h'
ax.set_title(title, fontsize=14,color='Black')
ax.plot(prob_qmapped_fullfield_dressed[0:-1:3,0:-1:3], \
    prob_qmapped_fullfield[0:-1:3,0:-1:3],'r.', markersize=0.2)    
ax.plot([0,1],[0,1],color='Black',lw=1.5)    
ax.set_ylim(0.0,1.0)
ax.set_xlim(0.0,1.0)
ax.set_xticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax.set_yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax.set_xlabel('Probability of > '+cthresh+' mm, quantile mapped, weighted & dressed', fontsize=10)
ax.set_ylabel('Probability of > '+cthresh+' mm, quantile mapped, weighted', fontsize=10)
ax.grid(color='LightGray',lw=0.5)

plot_title = 'probability_scatter_'+cthresh+'mm_'+\
    cyyyymmddhh+'_lead'+clead+'.png'
fig.savefig(plot_title, dpi=300)
plt.close()
print ('saving plot to file = ',plot_title)

    
