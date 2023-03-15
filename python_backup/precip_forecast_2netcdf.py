""" this script is designed to read in GEFSv12 forecast grib data, extract the CONUS subset, 
    and save them to a netcdf file.  Coded by Tom Hamill, Sep 2021, 
    tom.hamill@noaa.gov """
    
import numpy as np
import sys
import os
import time as timey
from netCDF4 import Dataset
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE


# =====================================================================

def set_domain_boundaries(cdomain):

    """ used grib file of 2.5-km blend output grid to determine bounding 
        lat and lon, and from that, the domain bounding indices for the 
        0.25 GEFSv12 reforecast data that will encompass the domain.    
    """
    if cdomain == 'conus': 
        jmin = 93
        jmax = 246
        imin = 368
        imax = 686
    elif cdomain == 'pr':
        jmin = 243
        jmax = 256
        imin = 649
        
        imax = 667   
    elif cdomain == 'ak':
        jmin = 19
        jmax = 161
        imin = 201
        imax = 967
    else:
        print ('invalid domain.  Exiting.')     
        sys.exit()    

    return jmin, jmax, imin, imax

# =====================================================================================

clead = sys.argv[1]  # beginning lead time in hours
cmonth = sys.argv[2] # Jan, Feb, etc.
cdomain = 'conus'
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)
master_directory = '/Volumes/NBM/'+cdomain+'_gefsv12/precip/netcdf/'

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time + surrounding months.  
#      All members, dates for this month have been
#      smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.

ncfile = master_directory + cmonth + '_apcp_sfc_h' + clead + '.nc'
print ('reading ',ncfile)
nc = Dataset(ncfile)
precip_forecast = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
nsamps, ny_gefsv12, nx_gefsv12 = np.shape(precip_forecast)
lons_1d = nc.variables['lons_fcst'][imin:imax]
lats_1d = nc.variables['lats_fcst'][jmin:jmax]
nlons = len(lons_1d)
nlats = len(lats_1d)
mem_num = nc.variables['mem_num'][:]
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
yyyymmddhh_fcst = nc.variables['yyyymmddhh_fcst'][:]
nc.close()
print ('done reading.')

# ---- set up for writing to output netCDf file

outfilename = master_directory + cmonth+'_'+cdomain+'_reforecast_precip_h'+clead+'.nc'
print ('writing netCDF data to ',outfilename)
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

xf = rootgrp.createDimension('xf',nlons)
xvf = rootgrp.createVariable('xf','f4',('xf',))
xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid"
xvf.units = "n/a"

yf = rootgrp.createDimension('yf',nlats)
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

mem_numv = rootgrp.createVariable('mem_num','i4',('sample',))
mem_numv.longname = "Member number (0-4)"

yyyymmddhh_init2 = rootgrp.createVariable('yyyymmddhh_init','i4',('sample',))
yyyymmddhh_init2.longname = "Initial condition date/time in yyyymmddhh format"

yyyymmddhh_fcst2 = rootgrp.createVariable('yyyymmddhh_fcst','i4',('sample',))
yyyymmddhh_fcst2.longname = \
    "Forecast date/time in yyyymmddhh format, end of 6-h precipitation accumulation period"

apcp_fcst = rootgrp.createVariable('apcp_fcst','f4',('sample','yf','xf',),
    zlib=True,least_significant_digit=6)
apcp_fcst.units = "mm"
apcp_fcst.long_name = "Ensemble precipitation forecast (mm)"
apcp_fcst.valid_range = [0.,10000.]
apcp_fcst.missing_value = np.array(-9999.99,dtype=np.float32)
        
# ---- metadata

rootgrp.stream = "s4" # ????
rootgrp.title = "GEFSv12 precipitation accumulation over the "+cdomain
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created Sep 2021 by Tom Hamill"
rootgrp.institution = \
    "Reforecast from ERSL/PSL using NCEP/EMC GEFSv12 circa 2020"
rootgrp.platform = "Model"
rootgrp.references = "https://noaa-gefs-retrospective.s3.amazonaws.com/index.html"

# ---- write variables

print ('writing data to output files.')
xvf[:] = range(nlons)
yvf[:] = range(nlats)
lonsf[:] = lons_1d[:]
latsf[:] = lats_1d[:]
yyyymmddhh_init2[:] = yyyymmddhh_init[:] 
yyyymmddhh_fcst2[:] = yyyymmddhh_fcst[:] 
samplev[:] = range(nsamps)
mem_numv[:] = mem_num[:]
for isamp in range(nsamps):
    apcp_fcst[isamp] = precip_forecast[isamp,:,:]

rootgrp.close()
print ('done!')

