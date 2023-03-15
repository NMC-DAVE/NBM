"""

python mswep_to_netCDF.py cmonth

For chosen month (01 to 12), extract 6-hourly files of MSWEP data, 
interpolate to CONUS NDFD grid and save to a new netCDF file. Do this
for all 6-hourly data for the chosen month and for years 2002-2019.

Dateutils is a Jeff Whitaker library for handling dates in yyyymmddhh format

Tom Hamill

"""

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

# ---- get the month and end time from the commmand line.  The first 00
#      hour analysis of the month will need to access the data from
#      the previous month.

cmonth = sys.argv[1] # 01 etc
imonth = int(cmonth) - 1

if cmonth == '01':
    cmonth_before = '12'
else:
    imonth_before = int(cmonth)-1
    if imonth_before < 10:
        cmonth_before = '0'+str(imonth_before)
    else:
        cmonth_before = str(imonth_before)
        
daysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
daysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
day1900 = datetime(1900,1,1,0)
hours1900 = dateto_hrs_since_day1CE(day1900, mixedcal=True)

# ---- get the lat/lons of the output NDFD CONUS grid by reading in
#      a sample file with the data.   These should be oriented S to N as 
#      the Basemap interp requires.

infile = 'ccpa.t00z.06h.ndgd2p5.conus.gb2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select()[0]
lats_ndfd, lons_ndfd = fcst.latlons()
if lats_ndfd[0,0] > lats_ndfd[-1,0]: 
    flipud = True
else:
    flipud = False
if flipud == True:
    lats_ndfd = np.flipud(lats_ndfd)
    lons_ndfd = np.flipud(lons_ndfd)
nlats_ndfd, nlons_ndfd = np.shape(lons_ndfd)
flatlon.close()
print ('min, max lons_ndfd = ', np.min(lons_ndfd), \
    np.max(lons_ndfd))
       
zeros = np.zeros((nlats_ndfd, nlons_ndfd), dtype=np.float32)
ones = np.ones((nlats_ndfd, nlons_ndfd), dtype=np.float32)

# --- read in the land-water mask

infile = 'ndfd_terrain_landwater.grib2.gb2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select()[0]
landmask = fcst.values
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

infile = 'various_nbm_plus_mask.nc'
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

# ---- read in a sample mswep lat/lon to see if need to flip
#      so that Basemap interp works.

mswep_output_directory = '/Volumes/Backup Plus/mswep/'
infile = mswep_directory + '200001_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
lats_mswep = nc.variables['lats'][:,:]
if lats_mswep[0,0] > lats_mswep[-1,0]: 
    flipud_mswep = True
else:
    flipud_mswep = False
nc.close()

mninenine = -99.99*np.ones((nlats_ndfd, nlons_ndfd), dtype=np.float64)

# ---- process all years for this month

for iyear in range(2002,2020):
    
    cyear = str(iyear)
    
    print ('****** processing year = ', iyear)
    
    # ---- determine the days of the month
    if iyear%4 == 0:
        ndays = daysomo_leap[imonth]
    else:
        ndays = daysomo[imonth]
    
    # ---- open netCDF output file and deal with all the variable definition and 
    #      such.
     
    outfile = '/Volumes/NBM/conus_panal/'+cyear+cmonth+'_mswep_on_ndfd_grid_6hourly.nc'
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
    
    conusmask = ncout.createVariable('conusmask','i4',('yf','xf',))
    latsa.long_name = "mask (1=land, 0=water)"
    latsa.units = "none"

    yyyymmddhh_begin = ncout.createVariable('yyyymmddhh_begin','i4',('time',))
    yyyymmddhh_begin.longname = \
        "Precip accumulation period beginning in yyyymmddhh format"

    yyyymmddhh_end = ncout.createVariable('yyyymmddhh_end','i4',('time',))
    yyyymmddhh_end.longname = \
        "Precip accumulation period ending in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    apcp_anal = ncout.createVariable('apcp_anal','f4',('time','yf','xf',),
        zlib=True,least_significant_digit=4)
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
    conusmask[:] = finalmask[:,:]

    # ---- metadata

    ncout.title = "NDFD CONUS domain, 6 hourly accum."
    ncout.history = " "
    ncout.institution =  "NCEP/EMC"
    ncout.platform = "Precipitation analysis"
    ncout.references = "DOI: 10.1175/JHM-D-11-0140.1"

    # ---- loop thru all dates, read reforecasts, and munge them into netCDF...

    ktr = 0
    for iday in range(1,ndays+1):
        
        infile2 = ''
        if iday < 10:
            cday = '0'+str(iday)
        else:
            cday = str(iday)
            
        cyyyymmdd = cyear + cmonth + cday
        for chour in ['00','06','12','18',]:
            
            precip_final = np.zeros((nlats_ndfd, nlons_ndfd), dtype=np.float64)
                
            ihour = int(chour)
            if ihour == 0:
                chour_begin = '18'
                chour_end = '00'
            elif ihour == 6:
                chour_begin = '00'
                chour_end = '06'
            elif ihour == 12:
                chour_begin = '06'
                chour_end = '12'
            elif ihour == 18:
                chour_begin = '12'
                chour_end = '18'
        
            cyyyymmddhh_end = cyear + cmonth + cday + chour
            cyyyymmddhh_begin = dateshift(cyyyymmddhh_end, -6)
            iyyyymmddhh = int(cyyyymmddhh_end)

            # ---- read the MSWEP data
                
            if iday == 1 and chour_end == '00':
                
                # the 18-00 UTC data on the first day of month should be stored with 
                # the data from the previous month.   For example, if
                # iyyyymmddhh == 2002010100, the data period is 2001123118 to
                # 2002010100 and it'd be stored in the Dec 2001 netCDF file
                
                if cmonth == '01':
                    cyear_before = str(int(cyear)-1)
                    infile2 = mswep_directory + cyear_before + \
                        cmonth_before + '_on_ndfd_grid_6hourly.nc'
                else:
                    infile2 = mswep_directory + cyear + \
                        cmonth_before + '_on_ndfd_grid_6hourly.nc'
            else:
                infile2 = mswep_directory + cyear + cmonth + '_on_ndfd_grid_6hourly.nc'
            nc = Dataset(infile2)
            print (infile2)
            yyyymmddhh_end_in = nc.variables['yyyymmddhh_end'][:]
            idx = int(np.where(yyyymmddhh_end_in == iyyyymmddhh)[0])
            precip_mswep = nc.variables['apcp_anal'][idx,:,:].astype('d')
            if flipud_mswep == True:
                precip_mswep = np.flipud(precip_mswep)
            nc.close()        
            
            # ---- write record to file.
            
            yyyymmddhh_begin[ktr] = int(cyyyymmddhh_begin)
            yyyymmddhh_end[ktr] = int(cyyyymmddhh_end)
            apcp_anal[ktr] = precip_mswep[:,:]
            ktr = ktr + 1

    ncout.close()
    print ('writing to ', outfile, ' completed.')






