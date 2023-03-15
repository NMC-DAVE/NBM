"""
plot_average_gain_2018.py

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
from read_reanalysis_timeseries import read_reanalysis_timeseries
from reformat_2d_to_4d_f90 import reformat_2d_to_4d_f90

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'
#rcParams['legend.fontsize']='large'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
cwarm = sys.argv[2]
ilead = int(clead)


if clead == '24':
    efold_random = '400.0'
    efold_bias = '1400.0'
elif clead == '120':
    if cwarm == 'warm':
        efold_random = '400.0'
        efold_bias = '1400.0'
    else:
        efold_random = '800.0' # '800.0'
        efold_bias = '1400.0' # '1400.0'
else:
    print ('invalid lead time. Stopping.')
    sys.exit()
        

cpath_Bx = '/Volumes/Backup Plus/ecmwf/Bx/'
cpath_Bbeta = '/Volumes/Backup Plus/ecmwf/Bbeta/'
cpath_gain = '/Volumes/Backup Plus/python/KFgain/'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
date_list_forecast = ['2018010100']

# --- read in a single reanalysis in order to get lats and lons.  Find the 
#     index of the nearest grid point.

analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
nlats, nlons = np.shape(lons)
print (nlats,nlons)
print ('min, max lons ', np.min(lons), np.max(lons))
print ('min, max lats ', np.min(lats), np.max(lats))
lons = lons - 360.
lats_1d = lats[:,0]
lons_1d = lons[0,:]

# --- read the Kalman filter gain for bias

gain_infile = cpath_gain + '2018_KFgain_flocal'+str(efold_random)+\
    '_blocal'+str(efold_bias)+'_2018_'+cwarm+'_lead'+\
    clead+'.cPick'
inf = open(gain_infile, 'rb')
Kalman_gain_beta_4d = cPickle.load(inf)
inf.close()
               
# ---- compute average of 1-point correlation maps

Kalman_gain_avg_map = np.zeros((41,81), dtype=np.float64)
ktr = 0
for ilocn in range(41,78,3):
    for jlocn in range(21,38,3):
        print (jlocn, ilocn)
        ktr = ktr+1
        Kalman_gain_beta_map = Kalman_gain_beta_4d[:,:,jlocn,ilocn]
        Kalman_gain_avg_map[:,:] = Kalman_gain_avg_map[:,:] + \
            Kalman_gain_beta_map[jlocn-20:jlocn+21,ilocn-40:ilocn+41]
                        
Kalman_gain_avg_map = Kalman_gain_avg_map / float(ktr)        
Kmax = np.max(np.abs(Kalman_gain_avg_map ))                      
Kalman_gain_avg_map = Kalman_gain_avg_map  / Kmax          
                             
# --- now contour the normalized Kalman average 

clevs = [-1.0,-0.7,-0.5,-0.3,-0.1,-0.05,0.05,0.1,0.3,0.5,0.7,1.0]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']    
    
fig = plt.figure(figsize=(6,3.8))
axloc = [0.12,0.27,0.85,0.62]
ax1 = fig.add_axes(axloc)
title = 'Average normalized Kalman gain map,\n'+cwarm+\
    '-season 2018, '+clead+'-hour forecast'
ax1.set_title(title, fontsize=11,color='Black')
ax1.set_xlabel('east-west distance (degrees)')
ax1.set_ylabel('north-south distance (degrees)')
ax1.set_xticks(range(-40,41,5))
ax1.set_yticks(range(-20,21,5))
ax1.grid(True, lw=0.25, color='Gray')


CS2 = ax1.contourf(range(-40,41), range(-20,21), \
    Kalman_gain_avg_map, clevs,\
    cmap=None, colors=colorst ,extend='both')
ax1.plot([0,0],[0,0],marker='.',markersize=5,color='Black')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
#cax = divider.append_axes("bottom", size="3%", pad=.4)
cax = fig.add_axes([0.05,0.11,0.9,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('average Kalman gain / max(absolute value(average Kalman gain))',fontsize=7)

# ---- set plot title

plot_title = 'average_normalized_Kalman_gain_'+cwarm+'2018_f'+\
    clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')

