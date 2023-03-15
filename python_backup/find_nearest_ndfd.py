import os, sys
import numpy as np
from netCDF4 import Dataset

clon = sys.argv[1]
clat = sys.argv[2]
rlon = float(clon)
rlat = float(clat)

master_directory_panal_spline = '/Volumes/NBM/conus_panal/CDF_spline/'

infile = master_directory_panal_spline +  \
    'Jan_conus_CCPA_spline_info_h00UTC.nc'
nc = Dataset(infile)
lons_ndfd = nc.variables['lons'][:,:]
print ('min, max lons_ndfd = ', np.min(lons_ndfd), np.max(lons_ndfd))
lats_ndfd = nc.variables['lats'][:,:]
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
nc.close() 

distx = np.abs(lons_ndfd - rlon)
disty = np.abs(lats_ndfd - rlat)
dist = distx**2 + disty**2

i = np.argmin(dist)
print (i)

jmin, imin = np.unravel_index(i, shape=(ny_ndfd, nx_ndfd), order='C')
print ('jmin, imin, lat, lon = ', jmin, imin, lats_ndfd[jmin, imin], lons_ndfd[jmin,imin])