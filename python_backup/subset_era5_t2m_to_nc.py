"""
subset_gefsv12_t2m_to_nc.py

for input lead time, develop a netCDF file of the time series of 
reforecasts for this lead time.

coded by: Tom Hamill, Mar 2021, tom.hamill@noaa.gov 

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from dateutils import daterange, dateshift, datetohrs
import pygrib # grib-reading routine
import scipy.stats as stats

# =====================================================================    

def get_domain_subset_of_era5(cyyyymmddhh, ifirst):

    """ read in the global forecast.  Subset to CONUS. """

    nib = 886  # precomputed for set of lat-lons that surrounds
    nie = nib+318 # the NDFD CONUS 2.5-km Lambert conformal domain.
    njb = 131 
    nje = njb+153 
    
    # ---- read in forecast grid covering the whole globe.
    
    input_directory = '/Volumes/Backup Plus/era5/t2m/'
    infile = input_directory + cyyyymmddhh + \
        '_ge'+cmem+'.t00z.pgrb2s.0p25.f' + clead 
        
    nc = Dataset(infile)


    t2m = np.squeeze(nc.variables['air'][idx,:,:])
    if ifirst == True:
        lons_in = nc.variables['lons'][:,:]
        lats_in = nc.variables['lats'][:,:]  
        time = nc.variables['time'][:]
        ntimes = len(times)
        yyyymmddhh = np.zeros((ntimes), dtype=np.int32)
        hourssince1CE_2000 = datetohrs('2000010100')
        for itime in enumerate(times):
            
            
            time1 = ncin.variables['time'][ifirst]
            hours1 = hours1900+time1*24 #hrs_since_day1CE_todate(hours1900+time1*24, mixedcal=True)
            #print (   'hours1 = ', hours1)
            cyyyymmddhh_begin = hrstodate(hours1,mixedcal=True)

            time2 = ncin.variables['time'][isecond]
            hours2 = hours1900+time2*24 + 3 #hrs_since_day1CE_todate(hours1900+time2*24 + 3, mixedcal=True)
            #print (   'hours2 = ', hours2)
            cyyyymmddhh_end = hrstodate(hours2,mixedcal=True)
            print ('   iday, ktr, cyyyymmddhh_begin, end = ', iday, ktr, \
                cyyyymmddhh_begin, cyyyymmddhh_end)
            
            
            datetime_list = []
            for i in range(len(dates)):
                dates[i] = int(strdates[i])
                hourssince2000[i] = datetohrs(strdates[i]) - hourssince1CE_2000
                datetime_list.append(datetime.strptime(strdates[i], "%Y%m%d%H"))
                decimalyear[i] = 2000. + hourssince2000[i]/(24.*365.25)
            timediff = datetime_list[len(dates)-1] - datetime_list[0]
            
            
              
        
    print ('   reading from ', infile)
    endStep = int(clead)
    istat, t2m_input, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    # ---- subset for CONUS.
    
    t2m = t2m_input[njb:nje,nib:nie]
    lons_1D = lons_full[0,nib:nie]-360.
    lats_1D = lats_full[njb:nje,0]

    return t2m, lons_1D, lats_1D
    
# =====================================================================    

def initialize_netCDF(master_directory_ncout, clead, nx, ny, lons_in, lats_in):

    """ initialize the output netCDF file """

    outfile = master_directory_ncout + 't2m_conus_'+clead+'h.nc'
    print ('initializing netCDF file ',outfile)
    nc_fullfield = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    
    # --- initialize dimensions, variable names
        
    xf = ncout.createDimension('xf',nx)
    xvf = ncout.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny)
    yvf = ncout.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number on 1/4-degree lat-lon grid"
    yvf.units = "n/a"

    sample = ncout.createDimension('sample',None)
    samplev = ncout.createVariable('samplev','i4',('sample',))
    samplev.units = "n/a"

    lons_out = ncout.createVariable('lons','f4',('yf','xf',))
    lons_out.long_name = "longitude"
    lons_out.units = "degrees_east"

    lats_out = ncout.createVariable('lats','f4',('yf','xf',))
    lats_out.long_name = "latitude"
    lats_out.units = "degrees_north"

    yyyymmddhh_init = ncout.createVariable('yyyymmddhh_init','i4',('sample',))
    yyyymmddhh_init.longname = "Initial condition date/time in yyyymmddhh format"

    yyyymmddhh_fcst = ncout.createVariable('yyyymmddhh_fcst','i4',('sample',))
    yyyymmddhh_fcst.longname = "Forecast valid date/time in yyyymmddhh format"
    
    memnum = ncout.createVariable('member','i4',('sample',))
    memnum.longname = "member number"

    # --- declare the single-level variable information on lat-lon grid

    t2m = ncout.createVariable('t2m','f4',\
        ('sample','yf','xf',),zlib=True,least_significant_digit=3)
    t2m.units = "n/a"
    t2m.long_name = '2-meter temperature in deg K.'
    t2m.valid_range = [220., 320.]
    t2m.missing_value = np.array(-99.99,dtype=np.float32)
    
    # ---- metadata

    ncout.title = 'GEFSv12 reforecast t2m forecast data for this lead on 0.25-deg'+\
        ' grid surrounding the CONUS'
    ncout.history = "GEFSv12 implemented at NCEP/EMC Sep 2020"
    ncout.institution =  "NCEP/EMC and PSL"
    ncout.platform = ""
    ncout.references = ""

    # ---- initialize

    xvf[:] = np.arange(nx)
    yvf[:] = np.arange(ny)
    lons_out[:] = lons_in[:,:]
    lats_out[:] = lats_in[:,:]
    
    return ncout

# =====================================================================    
    
def write_t2m_to_netcdf (isamp, imem, cyyyymmddhh, \
    cyyymmddhh_fcst, t2m_out):    
    
    """ write the t2m record to netCDF file """
    
    # ---- write probabilities and close

    yyyymmddhh_init[isamp] = int(cyyyymmddhh)
    yyyymmddhh_fcst[isamp] = int(cyyyymmddhh_fcst)
    t2m[isamp] = t2m[:,:]
    memnum[isamp] = imem

    istat = 0
    return istat
    
# =====================================================================    

clead = sys.argv(1)  # 3 digits


master_directory_ncout = '/Volumes/NBM/gefsv12/t2m/conus_netCDF/'
date_list = daterange('2000010100','2019123100',24)
isamp = 0
for idate, cyyyymmddhh in enumerate(date_list):
    cyyyymmddhh_fcst = dateshift(cyyyymmddhh, int(clead))
    for imem, cmem in enumerate(['c00','p01','p02','p03','p04']):
        
        # --- read in grib file t2m data for this date and lead time
        
        t2m, lons_1D, lats_1D =  get_domain_subset_of_gefs(date, cmem, clead)

        # --- first time through?  Initialize output netCDF file.
        
        if idate == 0 and imem == 0:
            nx = len(lons_1D)
            ny = len(lats_1D)
            ncout = initialize_netCDF(master_directory_ncout, clead, \
                nx, ny, lons_1D, lats_1D)
                
        # --- write this netCDF record
        
        istat = write_t2m_to_netcdf (isamp, imem, cyyyymmddhh, \
            cyyymmddhh_fcst, t2m_out)
            
        isamp = isamp + 1
    
ncout.close() 

        