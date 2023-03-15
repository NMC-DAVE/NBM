"""

python ccpa_to_netCDF.py cmonth

For chosen month (01 to 12), extract grib files of CCPA 
on CONUS NDFD grid and save to a new netCDF file.

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
    dateto_hrs_since_day1CE, hrstodate, datetohrs, dateshift
from mpl_toolkits.basemap import Basemap, interp

# ---- get the month and end time from the commmand line

cmonth = sys.argv[1] # 01 etc
imonth = int(cmonth) - 1
daysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
daysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
day1900 = datetime(1900,1,1,0)
hours1900 = dateto_hrs_since_day1CE(day1900, mixedcal=True)


# ---- get the lat/lons of the output NDFD CONUS grid.   These are
#      oriented S to N as interp requires

infile = '/Volumes/Backup Plus/ccpa/ccpa.20180101/00/ccpa.t00z.06h.ndgd2p5.conus.gb2'
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
print ('min, max lons_ndfd = ', np.min(lons_ndfd), np.max(lons_ndfd))
#lons_ndfd = lons_ndfd  # convert to deg E.    
       
zeros = np.zeros((nlats_ndfd, nlons_ndfd), dtype=np.float32)


# ---- read in the CONUS mask.  Not sure about accuracy.

infile = '/Volumes/Backup Plus/ccpa/supplemental_locations_ndfd2p5_Jan.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conusmask'][:,:]
nc.close()

# ---- process all years for this month

#for iyear in range(2002,2020):
#for iyear in range(2002,2020):
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
     
    outfile = '../ccpa/'+cyear+cmonth+'_ccpa_on_ndfd_grid_6hourly.nc'
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
    conusmask[:] = conusmask_in[:,:]

    # ---- metadata

    ncout.title = "NDFD CONUS domain interpolated from CCPA, 6 hourly accum."
    ncout.history = "Interpolated CCPA provided by Yan Luo, NCEP/EMC, Dec 2020"
    ncout.institution =  "NCEP/EMC"
    ncout.platform = "Precipitation analysis"
    ncout.references = "DOI: 10.1175/JHM-D-11-0140.1"

    # ---- loop thru all dates, read reforecasts, and munge them into netCDF...

    ktr = 0
    for iday in range(1,ndays+1):
        
        if iday < 10:
            cday = '0'+str(iday)
        else:
            cday = str(iday)
            
        cyyyymmdd = cyear + cmonth + cday
        for chour in ['00','06','12','18',]:
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
            print (cyyyymmddhh_end)
            cyyyymmddhh_begin = dateshift(cyyyymmddhh_end, -6)
            print (cyyyymmddhh_begin, cyyyymmddhh_end)
        
            # --- read the 6-hourly CCPA file. Make sure none subzero.
            
            try:
                infile = '/Volumes/Backup Plus/ccpa/ccpa.'+cyyyymmdd+\
                    '/'+chour+'/ccpa.t'+chour+'z.06h.ndgd2p5.conus.gb2'
                print (infile)
                grb = pygrib.open(infile)
                panal = grb.select()[0]
                precip_ccpa = panal.values
                grb.close()
                precip_ccpa = np.where(precip_ccpa < 0.0, zeros, precip_ccpa)
                if flipud == True:
                    precip_ccpa = np.flipud(precip_ccpa)
                    
                # ---- save to netCDF file.

                timev[ktr]             = ktr  
                yyyymmddhh_begin[ktr]  = int(cyyyymmddhh_begin)
                yyyymmddhh_end[ktr]    = int(cyyyymmddhh_end)
                apcp_anal[ktr]         = precip_ccpa
                ktr = ktr + 1
                
            except:
                print ('whoops!   some problem with ', infile)

    ncout.close()






