# examine_precipitation_timeseries.py cmonth cyyyymmddhh clon clat

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats


# =====================================================================

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx
    
# =====================================================================

cmonth = sys.argv[1]
cyyyymmddhh = sys.argv[2]
clon = sys.argv[3]
clat = sys.argv[4]
rlon_input = float(clon)
rlat_input = float(clat)

cleads = ['006','012','018','024','030',\
    '036','042','048','054','060',\
    '066','072','078','084','090',\
    '096','102','108','114','120',\
    '126','132','138','144','150',\
    '156','162','168','174','180',\
    '186','192','198','204','210',\
    '216','222','228','234','240']
nleads = len(cleads)
    
jmin = 93 # conus
jmax = 246
imin = 368
imax = 686
ny = jmax - jmin
nx = imax - imin
zeros = np.zeros((ny,nx),dtype=np.float64)


# ---- open the appropriate 240-h file and read the forecast for this case day

master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
cmonth = cmonths[imm]
ncfile = master_directory + cmonth+ \
    '_conus_reforecast_ens_mean_0_to_240h.nc'
nc = Dataset(ncfile)
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
idx = np.where(yyyymmddhh_init == int(cyyyymmddhh))[0]
print ('240-h idx = ', idx)
apcp_fcst = np.squeeze(nc.variables['apcp_fcst'][idx,:,:])
lons_1d = nc.variables['lons_fcst'][:]
lats_1d = nc.variables['lats_fcst'][:]
ny = len(lats_1d)
nx = len(lons_1d)
nc.close()


# ---- read in all the lead times for this month and get the 
master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
nsamps = 0
print ('lead (h)   member 1   member 2   member 3   member 4   member 5    mean ')
for clead in cleads:
    
    # --- read in all the data for this lead time, and set up arrays 
    #     first time through.
    
    ncfile = master_directory + cmonth + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    
    if clead == '006':
        
        yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
        idx = np.where(yyyymmddhh_init == int(cyyyymmddhh))[0]
        print (idx, yyyymmddhh_init[idx])
        ncases = len(yyyymmddhh_init) // 5
        lons_1d = nc.variables['lons_fcst'][imin:imax]
        lats_1d = nc.variables['lats_fcst'][jmin:jmax]
        if np.max(lons_1d) > 180.: lons_1d = lons_1d - 180.
        lons_2d, lats_2d = np.meshgrid(lons_1d, lats_1d)
        
        # ---- find index of nearest grid point to chosen lon/lat

        print ('rlon_input, rlat_input = ', rlon_input, rlat_input)
        distx = np.abs(lons_2d - rlon_input)
        disty = np.abs(lats_2d - rlat_input)
        dist = distx**2 + disty**2
        i = np.argmin(dist)
        print (i)
        jchoose, ichoose = np.unravel_index(i, shape=(ny, nx), order='C')    
        print (lons_2d[jchoose,ichoose], lats_2d[jchoose,ichoose])
        print ('jchoose, ichoose = ', jchoose, ichoose)
        
    precip_in = nc.variables['apcp_fcst'][idx,jchoose,ichoose]
    print (clead, precip_in, np.mean(precip_in), yyyymmddhh_init[idx[0]])
    nc.close()
    


