"""
python era5_renewables_convert_netcdf.py infile

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
from netCDF4 import Dataset
from dateutils import daterange, dateshift, datetohrs, hrstodate
import scipy.stats as stats
import datetime


# ------------------------------------------------------------------

def convert_time_to_yyyymmddhh_format(hours_since_1900):
    hours_1900_to_zero = datetohrs('1900010100')
    hours_sum = hours_1900_to_zero + hours_since_1900
    #hours_sum = [h+hours_1900_to_zero for h in hours_since_1900_list]
    yyyymmddhh = hrstodate(hours_sum,mixedcal=True)
    return yyyymmddhh

# ------------------------------------------------------------------


def write_netcdf_files(outfile, units, shortName, long_name, ny, nx, \
        latitude, longitude, data_to_output, ndates_check,\
        date_list_begin_thisyear, date_list_end_thisyear):
        
    """ initialize and write data to the output netCDF file """

    print ('  initializing netCDF file ',outfile)
    nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    # --- initialize dimensions, variable names

    xf = ncout.createDimension('xf',nx)
    xvf = ncout.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny)
    yvf = ncout.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number on 1/4-degree lat-lon grid"
    yvf.units = "n/a"

    sample = ncout.createDimension('sample',ndates_check)
    samplev = ncout.createVariable('samplev','i4',('sample',))
    samplev.units = "n/a"

    longitude_out = ncout.createVariable('longitude','f4',('xf',))
    longitude_out.long_name = "longitude"
    longitude_out.units = "degrees (negative is west)"

    latitude_out = ncout.createVariable('latitude','f4',('yf',))
    latitude_out.long_name = "latitude"
    latitude_out.units = "degrees north"

    yyyymmddhh_begin_out = ncout.createVariable('yyyymmddhh_begin','i4',('sample',))
    yyyymmddhh_begin_out.longname = "Beginning of analysis averaging period in yyyymmddhh format"

    yyyymmddhh_end_out = ncout.createVariable('yyyymmddhh_end','i4',('sample',))
    yyyymmddhh_end_out.longname = "End of analysis averaging period in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    output_variable = ncout.createVariable(shortName,'f4',\
        ('yf','xf','sample'),zlib=True,least_significant_digit=3)
    output_variable.units = units
    output_variable.long_name = long_name

    # ---- metadata

    ncout.title = 'ERA5 daily-averaged (00-21 UTC) '+shortName+' time series for this year'
    ncout.history = history
    ncout.institution =  "Copernicus, with data munging at PSL by Tom Hamill"

    # ---- initialize

    xvf[:] = np.arange(nx)
    yvf[:] = np.arange(ny)
    longitude_out[:] = longitude[:]
    latitude_out[:] = longitude[:]
    yyyymmddhh_begin_out[:] = int(date_list_begin_thisyear)
    yyyymmddhh_end_out[:] = int(date_list_end_thisyear)
    output_variable[:] = data_to_output[:,:,:]
    nc.close()
    istat = 0
    return istat

# ------------------------------------------------------------------

def read_netcdf_reanalysis_data_times_lat_lon(infile):
    
    nc = Dataset(infile)
    print ('reading from ', infile)
    times = nc.variables['time'][:]
    latitude = nc.variables['latitude'][:]
    latitude = nc.variables['longitude'][:]
    return times, latitude, longitude
    
# ------------------------------------------------------------------

# ---- read off the list of times

infile = sys.argv[1]
cyearb = infile[5:9]
cyeare = infile[10:14]
print ('cyearb, cyeare = ', cyearb, cyeare)

history = '2021-09-24 21:41:51 GMT by grib_to_netcdf-2.20.0: '+\
    '/opt/ecmwf/mars-client/bin/grib_to_netcdf -S param -o '+\
    '/cache/data6/adaptor.mars.internal-1632517596.466021-17844-12-7034697c-2681-43f0-9e13-295de7ef63b0.nc'+\
    '/cache/tmp/7034697c-2681-43f0-9e13-295de7ef63b0-adaptor.mars.internal-1632512130.8168957-17844-14-tmp.grib'
conventions = "CF-1.6"
infile = '/data/thamill/era5/'+infile
times, latitude, longitude = read_netcdf_reanalysis_data_times(infile):
nlats = len(latitude)
nlons = len(longitude)

# ---- convert the first and last times to yyyymmddhh format

yyyymmddhh_begin = convert_time_to_yyyymmddhh_format(times[0])
yyyymmddhh_end = convert_time_to_yyyymmddhh_format(times[-1])
print ('yyyymmddhh_begin, _end = ', yyyymmddhh_begin, yyyymmddhh_end)
yyyymmddhh_list = daterange(yyyymmddhh_begin, yyyymmddhh_end, 3)
ndates_3hourly = len(yyyymmddhh_list)
print ('expected dates = 35064.  Actual dates = ',ndates_3hourly)
print ('yyyymmddhh_list[0], [-1] = ', yyyymmddhh_list[0], yyyymmddhh_list[-1])
ndates_daily = ndates_3hourly // 8
print ('ndates_daily = ', ndates_daily)
u100_mean_dailies = np.zeros((ndates_daily, nlats, nlons), dtype=np.float32)
v100_mean_dailies = np.zeros((ndates_daily, nlats, nlons), dtype=np.float32)
ssrd_mean_dailies = np.zeros((ndates_daily, nlats, nlons), dtype=np.float32)
date_list_daily_begin = []
date_list_daily_end = []
print ('looping over dates in this netCDF file ...')
for idate, date in enumerate(yyyymmddhh_list[0:-1:8])

    print ('  processing idate, date = ', idate, date)
    
    # ---- read in the eight dates for this day
    
    u100 = nc.variables['u100'][idate:idate+8,:,:]
    v100 = nc.variables['v100'][idate:idate+8,:,:]
    ssrd = nc.variables['ssrd'][idate:idate+8,:,:]
            
    # ---- get the mean across times.
    
    u100_daily_mean = np.mean(u100, axis=0)
    v100_daily_mean = np.mean(v100, axis=0)
    ssrd_daily_mean = np.mean(ssrd, axis=0)
        
    # ---- plunk into output array
    
    u100_mean_dailies[idate,:,:] = u100_daily_mean[:,:]
    v100_mean_dailies[idate,:,:] = v100_daily_mean[:,:]
    ssrd_mean_dailies[idate,:,:] = ssrd_daily_mean[:,:]
    
    # ---- append date for beginning and end of averaging period to lists
    #      of dates
    
    date_list_daily_begin.append(date)
    date_list_daily_end.append(dateshift(date, 21))
    
nc.close()

# ---- setup for netCDF file output

print ('------ beginning 2nd step of writing out each year ------ ')
year_list = range(int(cyearb), int(cyeare)+1)
print ('writing outputs for years ', year_list)

for year in year_list:
    
    print ('processing year = ', year)    
    cyear = str(year)
    date_begin = cyear+'010100'
    date_end = cyear+'123100'    
    
    idx_begin = np.where(date_list_daily_begin == int(date_begin))[0]
    idx_end = np.where(date_list_daily_begin == int(date_end))[0]
    print ('   indices of beginning, ending dates for ',cyear,' = ',idx_begin, idx_end)
    date_list_begin_thisyear = daterange(date_begin, date_end, 24)
    date_list_end_thisyear = daterange(dateshift(date_begin, 21), dateshift(date_end,21), 24)
    print ('   date_list_begin_thisyear = ', date_list_begin_thisyear)
    print ('   date_list_end_thisyear = ', date_list_end_thisyear)
    
    # ---- extract this year's data
    
    u100_daily_mean_thisyear = u100_mean_dailies[idx_begin:idx_end+1,:,:]
    v100_daily_mean_thisyear = v100_mean_dailies[idx_begin:idx_end+1,:,:]
    ssrd_daily_mean_thisyear = ssrd_mean_dailies[idx_begin:idx_end+1,:,:]
    
    # ---- transpose the array to put the list of dates as the trailing dimension.
    #      in this way, for a given [i,j], all training samples across dates 
    #      will be sequential on disk.
    
    u100_daily_mean_thisyear_t = np.transpose(u100_daily_mean_thisyear, (1, 2, 0))
    v100_daily_mean_thisyear_t = np.transpose(u100_daily_mean_thisyear, (1, 2, 0))
    ssrd_daily_mean_thisyear_t = np.transpose(ssrd_daily_mean_thisyear, (1, 2, 0))
    nycheck, nxcheck, ndates_check = np.shape(u100_daily_mean_thisyear_t)

    # ----- write out the netCDF records for u100
    
    outfile = '/data/thamill/era5/u100_dailymean_'+cyear+'.nc'
    units = "m s**-1" 
    shortName = 'u100'
    long_name = "100 metre U wind component (daily average)" 
    istat = write_netcdf_files(outfile, units, shortName, long_name, ny, nx, \
        latitude, longitude, u100_daily_mean_thisyear_t, ndates_check,\
        date_list_begin_thisyear, date_list_end_thisyear)

    # ----- write out the netCDF records for v100
    
    outfile = '/data/thamill/era5/v100_dailymean_'+cyear+'.nc'
    units = "m s**-1" 
    shortName = 'v100'
    long_name = "100 metre V wind component (daily average)" 
    istat = write_netcdf_files(outfile, units, shortName, long_name, ny, nx, \
        latitude, longitude, v100_daily_mean_thisyear_t, ndates_check,\
        date_list_begin_thisyear, date_list_end_thisyear)
        
    # ----- write out the netCDF records for ssrd
    
    outfile = '/data/thamill/era5/ssrd_dailymean_'+cyear+'.nc'
    units = "J m**-2"
    shortName = 'ssrd'
    long_name = "Downward shortwave radiation at surface (daily average)" 
    istat = write_netcdf_files(outfile, units, shortName, long_name, ny, nx, \
        latitude, longitude, v100_daily_mean_thisyear_t, ndates_check,\
        date_list_begin_thisyear, date_list_end_thisyear)

