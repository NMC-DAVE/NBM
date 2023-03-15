"""

mswep_to_ndfd_grid.py cmonth

For chosen month, sum to 6-hourly and interpolate
MSWEP to CONUS grid and save to a new netCDF file.

Tom Hamill, Dec 2020

"""

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats
import pygrib
from netCDF4 import Dataset
from dateutils import hrs_since_day1CE_todate, \
    dateto_hrs_since_day1CE, hrstodate, datetohrs
from mpl_toolkits.basemap import Basemap, interp

# ---- get the month and end time from the commmand line

cmonth = sys.argv[1] # 01 etc
imonth = int(cmonth) - 1
daysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
daysomo_leap = [31, 29, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
day1900 = datetime(1900,1,1,0)
hours1900 = dateto_hrs_since_day1CE(day1900, mixedcal=True)


# ---- get the lat/lons of the output NDFD CONUS grid.   These are
#      oriented S to N as interp requires

infile = '/Volumes/Backup Plus/blend_domains/blend.t00z.qmd.f001.co.grib2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select()[0]
lats_ndfd, lons_ndfd = fcst.latlons()
nlats_ndfd, nlons_ndfd = np.shape(lons_ndfd)
flatlon.close()
print ('min, max lons_ndfd = ', np.min(lons_ndfd), np.max(lons_ndfd))
lons_ndfd = lons_ndfd + 360.  # convert to deg E.    
       
# ---- process all years for this month

for iyear in range(2000,2020):
#for iyear in range(2000,2001):
    cyear = str(iyear)
    
    print ('****** processing year = ', iyear)
    # ---- determine the days of the month
    if iyear%4 == 0:
        ndays = daysomo_leap[imonth]
    else:
        ndays = daysomo[imonth]
    
    # ---- read a test date and convert to current year/month/day/hour

    infile = '/Volumes/Backup Plus/mswep/'+cyear+cmonth+'.nc'
    print ('  reading from ',infile)
    ncin = Dataset(infile)
    #time = nc.variables['time'][:]  # time, after conversion, is start of 3-h period
    if iyear == 2000:
        lons_mswep = ncin.variables['lon'][:]
        print ('lons_mswep = ', lons_mswep)
        nlons_mswep = len(lons_mswep)
        lats_mswep = ncin.variables['lat'][:]
        nlats_mswep = len(lats_mswep)
        zeros = np.zeros((nlats_mswep, nlons_mswep), dtype=np.float32)
        if lats_mswep[0] > lats_mswep[-1]: # not oriented S to N so flip
            flipud = True
            lats_mswep = np.flipud(lats_mswep)
        lons_mswep = lons_mswep + 360. 
        #print ('lons_mswep = ', lons_mswep)   
        #sys.exit()

    # ---- open netCDF output file and deal with all the variable definition and 
    #      such.
     
    outfile = '../mswep/'+cyear+cmonth+'_on_ndfd_grid_6hourly.nc'
    print ('   writing to ',outfile)
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    xf = ncout.createDimension('xf',nlons_ndfd)
    xvf = ncout.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on NDFD grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',nlats_ndfd)
    yvf = ncout.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number on 1/4-degree lat-lon grid"
    yvf.units = "n/a"

    time = ncout.createDimension('time',None)
    timev = ncout.createVariable('time','f4',('time',))
    timev.units = "index to time dimension, that's all"

    lonsa = ncout.createVariable('lons','f4',('yf','xf',))
    lonsa.long_name = "longitude"
    lonsa.units = "degrees_east"

    latsa = ncout.createVariable('lats','f4',('yf','xf',))
    latsa.long_name = "latitude"
    latsa.units = "degrees_north"

    yyyymmddhh_begin = ncout.createVariable('yyyymmddhh_begin','i4',('time',))
    yyyymmddhh_begin.longname = \
        "Precip accumulation period beginning in yyyymmddhh format"

    yyyymmddhh_end = ncout.createVariable('yyyymmddhh_end','i4',('time',))
    yyyymmddhh_end.longname = \
        "Precip accumulation period ending in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    apcp_anal = ncout.createVariable('apcp_anal','f4',('time','yf','xf',),
        zlib=True,least_significant_digit=2)
    apcp_anal.units = "mm"
    apcp_anal.long_name = \
        "Interpolated 6-h accumulated MSWEP analysis on CONUS NDFD grid"
    apcp_anal.valid_range = [0.,1000.]
    apcp_anal.missing_value = np.array(-9999.99,dtype=np.float32)

    # ---- initialize

    xvf[:] = np.arange(nlons_ndfd)
    yvf[:] = np.arange(nlats_ndfd)
    lonsa[:] = lons_ndfd[:,:]
    latsa[:] = lats_ndfd[:,:]

    # ---- metadata

    ncout.title = "NDFD CONUS domain interpolated from MSWEP, 6 hourly accum."
    ncout.history = "Interpolated from MSWEP_V260 0.1 deg data Dec 2020 by Tom Hamill"
    ncout.institution =  "gloh2o.org"
    ncout.platform = "Precipitation analysis"
    ncout.references = "http://gloh2o.org/"

    # ---- loop thru all dates, read reforecasts, and munge them into netCDF...

    ktr = 0
    for iday in range(ndays):
        
        for chour in ['06','12','18','00']:
            ihour = int(chour)
            if ihour == 0:
                ifirst = 6 + iday*8
                chour_begin = '18'
                chour_end = '00'
            elif ihour == 6:
                ifirst = 0 + iday*8
                chour_begin = '00'
                chour_end = '06'
            elif ihour == 12:
                ifirst = 2 + iday*8
                chour_begin = '06'
                chour_end = '12'
            elif ihour == 18:
                ifirst = 4 + iday*8
                chour_begin = '12'
                chour_end = '18'
            isecond = ifirst + 1
            print ('ifirst, isecond = ', ifirst, isecond)
        
            # --- set the bounding times for the precipitation analyses, begin&end
        
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
        
            # --- read the two 3-hourly netCDF slices. Make sure none subzero.
         
            precip_first = ncin.variables['precipitation'][ifirst,:,:]
            precip_second = ncin.variables['precipitation'][isecond,:,:]
            precip_mswep = precip_first + precip_second
            precip_mswep = np.where(precip_mswep < 0.0, zeros, precip_mswep)
            if flipud == True:
                precip_mswep = np.flipud(precip_mswep)
        
            # ---- bilinear interpolate to the NDFD grid
        
            precip_ndfd = interp(precip_mswep, lons_mswep, lats_mswep, \
                lons_ndfd, lats_ndfd, checkbounds=False, masked=False, order=3)
                    
            # ---- save to netCDF file.


            yyyymmddhh_begin[ktr]  = int(cyyyymmddhh_begin)
            yyyymmddhh_end[ktr]    = int(cyyyymmddhh_end)
            apcp_anal[ktr]         = precip_ndfd
            ktr = ktr + 1

    ncin.close()
    ncout.close()






