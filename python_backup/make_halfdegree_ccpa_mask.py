# make_halfdegree_ccpa_mask.py

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats
import pygrib
from netCDF4 import Dataset
from dateutils import hrs_since_day1CE_todate, \
    dateto_hrs_since_day1CE, hrstodate, datetohrs, dateshift
from mpl_toolkits.basemap import Basemap, interp

# --- read in the land-water mask

infile = 'ndfd_terrain_landwater.grib2.gb2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select()[0]
landmask = fcst.values
ny, nx = np.shape(landmask)
ones = np.ones((ny,nx), dtype=int)
zeros = np.zeros((ny,nx), dtype=int)
landmask = np.where(landmask > 1.0, ones, zeros)
lats_landmask, lons_landmask = fcst.latlons()
if lats_landmask[0,0] > lats_landmask[-1,0]: 
    flipud = True
else:
    flipud = False
if flipud == True:
    lats_landmask = np.flipud(lats_landmask)
    lons_landmask = np.flipud(lons_landmask)
nlats_landmask, nlons_landmask = np.shape(lons_landmask)
flatlon.close()

# --- read in the valid CCPA mask.  This contains
#     Gulf and Atlantic points we don't want to use.
#     These will be filtered out later.
#     Got these masks from Eric Engle, MDL

infile = '/Volumes/Backup Plus/ccpa/various_nbm_plus_mask.nc'
nc = Dataset(infile)
lats_ccpa = nc.variables['latitude'][:,:]
lons_ccpa = nc.variables['longitude'][:,:]
validmask_ccpa = nc.variables['validmask_ccpa'][:,:]
if lats_ccpa[0,0] > lats_ccpa[-1,0]: 
    validmask_ccpa = np.flipud(validmask_ccpa)
nc.close()

# ---- make the final mask as landmask*validmask_ccpa. 
#      Through this, we will default to using the 
#      alternative MSWEP analysis at all water points.

finalmask = validmask_ccpa*landmask


lons_fcst = np.arange(-128.5, -57.51, 0.5)
lats_fcst = np.arange(57.5, 18.99, -0.5)
print (lats_fcst)
print (lons_fcst)