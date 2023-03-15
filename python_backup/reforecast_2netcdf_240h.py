# reforecast_2netcdf_240h.py

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats


cmonth = sys.argv[1]
cleads = ['006','012','018','024','030',\
    '036','042','048','054','060',\
    '066','072','078','084','090',\
    '096','102','108','114','120',\
    '126','132','138','144','150',\
    '156','162','168','174','180',\
    '186','192','198','204','210',\
    '216','222','228','234','240']
nleads = len(cleads)
    
jmin = 93
jmax = 246
imin = 368
imax = 686
ny = jmax - jmin
nx = imax - imin
zeros = np.zeros((ny,nx),dtype=np.float64)
    
# ---- read in all the lead times for this month and get the 
master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
nsamps = 0
for clead in cleads:
    
    # --- read in all the data for this lead time, and set up arrays 
    #     first time through.
    
    ncfile = master_directory + cmonth + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_in = nc.variables['apcp_fcst'][:,:,:]
    print ('np.shape(precip_in) = ', np.shape(precip_in) )
    if clead == '006':
        yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
        ncases = len(yyyymmddhh_init) // 5
        lons_1d = nc.variables['lons_fcst'][imin:imax]
        lats_1d = nc.variables['lats_fcst'][jmin:jmax]
        precip_sum = np.zeros((ncases,ny,nx), dtype=np.float64)
    nc.close()
    
    # ---- for each case day, sum over 5 reforecast members.
    
    for icase in range(ncases):
        if icase%20 == 0: print ('icase = ',icase)
        istart = icase*5
        #for imem in range(5):
        #    precip_sum[icase,:,:] = precip_sum[icase,:,:] + \
        #        np.sum(precip_in[istart:istart+5,jmin:jmax,imin:imax], axis=0)
        precip_sum[icase,:,:] = precip_sum[icase,:,:] + \
            np.sum(precip_in[istart:istart+5,jmin:jmax,imin:imax], axis=0)

# ---- after summing over all lead times, now divide the total sum by 
#      the number of reforecast members (5) to get the ensemble mean.

print('computing ensemble mean')
yyyymmddhh_init_out = np.zeros((ncases), dtype=np.int32)
for icase in range(ncases):
    precip_sum[icase,:,:] = precip_sum[icase,:,:] / 5.
    yyyymmddhh_init_out[icase] = yyyymmddhh_init[icase*5]
    #print(precip_sum[icase,61,58], yyyymmddhh_init_out[icase])

# ---- set up for writing to output netCDf file

outfilename = master_directory + cmonth+'_conus_reforecast_ens_mean_0_to_240h.nc'
print ('writing netCDF data to ',outfilename)
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

xf = rootgrp.createDimension('xf',nx)
xvf = rootgrp.createVariable('xf','f4',('xf',))
xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid"
xvf.units = "n/a"

yf = rootgrp.createDimension('yf',ny)
yvf = rootgrp.createVariable('yf','f4',('yf',))
yvf.long_name = "southward grid point number on 1/4-degree lat-lon grid"
yvf.units = "n/a"

lonsf = rootgrp.createVariable('lons_fcst','f4',('xf',))
lonsf.long_name = "longitude"
lonsf.units = "degrees_east"

latsf = rootgrp.createVariable('lats_fcst','f4',('yf',))
latsf.long_name = "latitude"
latsf.units = "degrees_north"

sample = rootgrp.createDimension('sample',None)
samplev = rootgrp.createVariable('samplev','i4',('sample',))
samplev.units = "n/a"

yyyymmddhh_init2 = rootgrp.createVariable('yyyymmddhh_init','i4',('sample',))
yyyymmddhh_init2.longname = "Initial condition date/time in yyyymmddhh format"

apcp_fcst = rootgrp.createVariable('apcp_fcst','f4',('sample','yf','xf',),
    zlib=True,least_significant_digit=6)
apcp_fcst.units = "mm"
apcp_fcst.long_name = "Ensemble-mean total 240-h precipitation forecast"
apcp_fcst.valid_range = [0.,10000.]
apcp_fcst.missing_value = np.array(-9999.99,dtype=np.float32)

# ---- metadata

rootgrp.stream = "s4" # ????
rootgrp.title = "GEFSv12 precipitation accumulation over the "+\
    "expanded MDL domain incl Guam, HI, AK, PR, CONUS"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created November 2020 by Tom Hamill"
rootgrp.institution = \
    "Reforecast from ERSL/PSL using NCEP/EMC GEFSv12 circa 2020"
rootgrp.platform = "Model"
rootgrp.references = "https://noaa-gefs-retrospective.s3.amazonaws.com/index.html"

# ---- initialize

xvf[:] = np.arange(nx)
yvf[:] = np.arange(ny)
lonsf[:] = lons_1d[:]
latsf[:] = lats_1d[:]

# ---- write netCDF records

for icase in range(ncases):
    apcp_fcst[icase] = precip_sum[icase,:,:] 
    yyyymmddhh_init2[icase] = yyyymmddhh_init_out[icase]
    
rootgrp.close()
print ('done!')
