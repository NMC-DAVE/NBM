""" this script will take yearly ERA5 data for a variable of choice, 
subset it to the CONUS domain, reorder the indices of the data, and 
write revised records back to a new netCDF file
"""



import numpy as np
import sys
import pygrib
import os
import time as timey
from netCDF4 import Dataset

from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp

def initialize_netCDF_metadata(variable)


    people = {'air': {'name': 'John', 'age': '27', 'sex': 'Male'},
              '': {'name': 'Marie', 'age': '22', 'sex': 'Female'}}


def initialize_output_netcdf_file(ny, nx, ndates, outfilename, \
    yyyymmddhh, variable_units, variable_name, variable_long_name, \
    valid_name, missing_value, least_significant_digit, title, \
    conventions, history, ):
    
    # ---- set up netCDF output file particulars

    print (outfilename)
    rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

    xf = rootgrp.createDimension('xf',ni)
    xvf = rootgrp.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid" 
    xvf.units = "n/a" 

    yf = rootgrp.createDimension('yf',None)
    yvf = rootgrp.createVariable('yf','f4',('yf',))
    yvf.long_name = "southward grid point number on 1/4-degree lat-lon grid" 
    yvf.units = "n/a" 

    time = rootgrp.createDimension('time',ndates)
    timev = rootgrp.createVariable('time','f4',('time',))
    timev.units = "hours since 1-1-1 00:00:0.0" 

    lons = rootgrp.createVariable('lons','f4',('xf',))
    lons.long_name = "longitude" 
    lons.units = "degrees (negative = west)" 

    lats = rootgrp.createVariable('lats','f4',('yf',))
    lats.long_name = "latitude" 
    lats.units = "degrees_north"  

    yyyymmddhh = rootgrp.createVariable('yyyymmddhh','i4',('time',))
    yyyymmddhh.longname = "Analysis date/time in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    analysis = rootgrp.createVariable(variable_name,'f4',('yf','xf','time'),
        zlib=True,least_significant_digit=least_significant_digit)  
    analysis.units = variable_units
    analysis.long_name = variable_long_name
    analysis.valid_range = valid_range
    analysis.missing_value = np.array(missing_value,dtype=np.float32)
    
    # ---- metadata 

    rootgrp.stream = ""
    rootgrp.title = title
    rootgrp.Conventions = conventions  
    rootgrp.history = history
    rootgrp.institution = "ECMWF/Copernicus"
    rootgrp.platform = "Model" 
    rootgrp.references = "https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3803"


# =====================================================================================


# ---- define bounding coordinates of the CONUS domain

nib = 
nie =
njb =
nje = 

ny = nje - njb 
nx = nie - nib
print ('ny, nx = ', ny, ny)





# ---- initialize

xvf[:] = np.arange(ni)
yvf[:] = np.arange(nj)

lonsf_out = \
    np.hstack((lons_quarterdeg_1d[nib1:nie1], lons_quarterdeg_1d[nib2:nie2]))
lonsf[:] = lonsf_out
print ('output lons from ',lonsf_out[0], lonsf_out[-1])
print ('output lats from ',lats_quarterdeg_1d[njb], lats_quarterdeg_1d[nje-1])
latsf[:] = lats_quarterdeg_1d[njb:nje] 
ensv[:] = range(5)



# ---- loop thru all dates, read reforecasts, and munge them into netCDF...

ktr = 0
work_array = np.zeros((nmembers,nj,ni),dtype=np.float)
zeros = np.zeros((nmembers,nj,ni),dtype=np.float)
work_array1 = np.zeros((nj,ni),dtype=np.float)
work_array2 = np.zeros((nj,ni),dtype=np.float)
variable = 'apcp_sfc_'
#print ('date_list = ',date_list)
ifirst = True
for idate, date_IC in zip(range(ndates), date_list):

    date_FC = dateshift(date_IC, ilead)
    print ('processing date_IC, date_FC = ', date_IC, date_FC)
    cyear = date_IC[0:4]
    chour = date_IC[8:10]
    cyyyymmddhh = date_IC
    cyearmo = date_IC[0:6]
    validityDate = int(date_FC)//100
    validityTime = (int(date_FC) - validityDate*100)*100
    #print ('   ',idate, 'initial = ',date_IC,' forecast = ',date_FC,\
    #    ' validityDate=',validityDate,' validityTime = ',\
    #    validityTime, timey.asctime())
    
    work_array[:,:,:] = 0.0
    work_array2[:,:] = 0.0
    work_array1[:,:] = 0.0
    
    if validityTime == 600 or validityTime == 1200 or \
        validityTime == 1800 or validityTime == 0:
        
        # --- to get 3-hourly precipitation accumulation, we'll 
        #     need to read the file 3 hours beforehand
        #     and subtract that amount
        
        for imem, cmem in enumerate(cmembers): 
            
            #print ('   processing member ', imem, ' ',cmem)
            filename, filename_d = \
                form_filename(cyear, date_IC, cmem, variable)
            istat2, grb_late = read_gribdata(filename_d, \
                filename, validityDate, validityTime)
            #print ('   read grib_late')
            
            date_FC_early = dateshift(date_FC, -3)
            validityDate_early = int(date_FC_early)//100
            validityTime_early = (int(date_FC_early) - \
                validityDate_early*100)*100
            istat2, grb_early = read_gribdata(filename_d, \
                filename, validityDate_early, validityTime_early)
            #print ('   read grib_early')
            
            grb_global = grb_late.values - grb_early.values
            work_array2[:,:] = np.hstack((grb_late.values\
                [njb:nje,nib1:nie1], grb_late.values[njb:nje,nib2:nie2]))
            work_array1[:,:] = np.hstack((grb_early.values\
                [njb:nje,nib1:nie1], grb_early.values[njb:nje,nib2:nie2]))
            work_array[imem,:,:] = work_array2[:,:] - work_array1[:,:]
            #print ('min, max work_array2 = ', np.min(work_array2), np.max(work_array2))
            #print ('min, max work_array1 = ', np.min(work_array1), np.max(work_array1))
            #print ('min, max work_array = ', np.min(work_array), np.max(work_array))
            #print ('   populated work arrays') 
    else: 
            
        # --- only need to read one file, no subtraction
        
        for imem, cmem in enumerate(cmembers): 
            #print ('   processing member ', imem, ' ',cmem)
            filename, filename_d = \
                form_filename(cyear, date_IC,  cmem, variable)
            istat, grb = read_gribdata(filename_d, \
                filename, validityDate, validityTime)
            #print ('   istat, grib files read')
            grb_global = grb.values
            work_array[imem,:,:] = np.hstack((grb_global\
                [njb:nje,nib1:nie1],grb_global[njb:nje,nib2:nie2]))
            #print ('   work array formed.')
            #print ('min, max work_array = ', np.min(work_array), np.max(work_array))
                 
    work_array = np.where(work_array < 0.0, zeros, work_array)
    apcp_fcst[ktr]  =  work_array[:,:,:]
    
    # ---- indicate the time since 0 AD, init and final time in yyyymmddhh format

    timev[ktr]             = datetohrs(date_FC)  # in our file, since 0 AD
    yyyymmddhh_init2[ktr]  = date_IC
    yyyymmddhh_fcst2[ktr]  = date_FC
    ktr = ktr + 1 

print ('writing to ',outfilename)
rootgrp.close()
istat = 0
return istat

