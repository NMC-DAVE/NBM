""" this script, for a list of dates and lead times, will 
extract the GEFS v12 3-hourly accumulated 
precip reforecast data and put it into a netCDF file 
"""

def reforecast_2netcdf (date_list, ilead, outfilename):

    import numpy as np
    import sys
    import pygrib
    import os
    import time as timey
    from netCDF4 import Dataset

    from dateutils import hrstodate, daterange, dayofyear, \
         splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
    from mpl_toolkits.basemap import interp


    # --- read grib data on a single level

    def read_gribdata(idxfilename, gribfilename, validityDate, validityTime):
        #grb=[]
        #print (idxfilename, gribfilename)
        istat = -1
        fexist_grib = False
        fexist_grib = os.path.exists(gribfilename)
        #fexist_idx  = False
        #fexist_idx  = os.path.exists(idxfilename)
        #print ('fexist_grib  = ', fexist_grib)
        #if fexist_grib and fexist_idx:
        if fexist_grib:
            try:
                fcstfile = pygrib.open(gribfilename)
                #fcstfile = pygrib.index(idxfilename)
                #print ('   q:',validityDate, validityTime)
                grb = fcstfile.select(shortName='tp',\
                    validityDate=validityDate, \
                    validityTime=validityTime)[0]
                gv = grb.values
                istat = 0
                fcstfile.close()
                return istat,grb
            except IOError:
                print ('   IOError in read_gribdata reading ', \
                    gribfilename, validityDate, validityTime)
                istat = -1
            except ValueError:
                print ('   ValueError in read_gribdata reading ', \
                    gribfilename, validityDate, validityTime)
                istat = -1
            except RuntimeError:
                print ('   RuntimeError in read_gribdata reading ', \
                    gribfilename, validityDate, validityTime)
                istat = -1
        return istat, grb


    # --- form the grib and grib index file name

    def form_filename(cyear, date_IC,  cmem, variable):
        master_directory = '/Volumes/Backup Plus/gefsv12/precip/'
        filename = master_directory + variable+date_IC+'_'+cmem+'.grib2'
        filename_idx = filename+'.idx'
        #print ('   ',filename)
        #print ('   ',filename_idx)
        istat = 0
        return filename, filename_idx

    # =====================================================================================

    # ---- read in a sample forecast in order to define the lat/lon information 
    #      for ~1/4-degree grid

    infile = '/Volumes/Backup Plus/gefsv12/precip/apcp_sfc_2019011000_p01.grib2'
    print (infile)
    flatlon = pygrib.open(infile)
    fcst = flatlon.select(shortName='tp', dataDate=20190110, validityTime=600)[0]
    lats_quarterdeg, lons_quarterdeg = fcst.latlons()
    nlats_quarterdeg, nlons_quarterdeg = lons_quarterdeg.shape
    lats_quarterdeg_1d = np.zeros((nlats_quarterdeg),dtype=np.float32)
    lons_quarterdeg_1d = np.zeros((nlons_quarterdeg),dtype=np.float32)
    lats_quarterdeg_1d[:] = lats_quarterdeg[:,0]
    lons_quarterdeg_1d[:] = lons_quarterdeg[0,:]
    lons_quarterdeg_1d = np.where(lons_quarterdeg_1d > 90.0, \
        lons_quarterdeg_1d - 360., lons_quarterdeg_1d)
    lats00 = lats_quarterdeg_1d[0]
    flatlon.close()

    #for ilon in range(nlons_quarterdeg):
    #   print (ilon, lons_quarterdeg_1d[ilon])
    #for ilat in range(nlats_quarterdeg):
    #    print (ilat, lats_quarterdeg_1d[ilat])
 
    # ---- define bounding coordinates on Gaussian grid for MDL domain
   
    nib1 = 518 # ~lon 220E
    nie1 = 1440 # up to ~lon 310E
    nib2 = 0 # ~lon 220E
    nie2 = 45 # up to ~lon 310E
    njb = 38 # lat ~ 80.5
    nje = 483 # down to lat ~ -30.5
    nj = nje - njb 
    ni1 = nie1 - nib1
    ni2 = nie2 - nib2
    ni = ni1+ni2 
    print ('nj, ni = ', nj, ni)

    # ---- initialize

    ndates = len(date_list)
    cmembers = ['c00','p01','p02','p03','p04']
    nmembers = len(cmembers)

    # ---- set up netCDF output file particulars

    print (outfilename)
    rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

    xf = rootgrp.createDimension('xf',ni)
    xvf = rootgrp.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid" 
    xvf.units = "n/a" 

    yf = rootgrp.createDimension('yf',nj)
    yvf = rootgrp.createVariable('yf','f4',('yf',))
    yvf.long_name = "southward grid point number on 1/4-degree lat-lon grid" 
    yvf.units = "n/a" 

    time = rootgrp.createDimension('time',None)
    timev = rootgrp.createVariable('time','f4',('time',))
    timev.units = "hours since 1-1-1 00:00:0.0" 

    ens = rootgrp.createDimension('ens',5)
    ensv = rootgrp.createVariable('ensv','i4',('ens',))
    ensv.long_name = "Ensemble member number (control, perts 1-4)" 
    ensv.units = " " 

    lonsf = rootgrp.createVariable('lons_fcst','f4',('xf',))
    lonsf.long_name = "longitude" 
    lonsf.units = "degrees_east" 

    latsf = rootgrp.createVariable('lats_fcst','f4',('yf',))
    latsf.long_name = "latitude" 
    latsf.units = "degrees_north"  

    yyyymmddhh_init2 = rootgrp.createVariable('yyyymmddhh_init','i4',('time',))
    yyyymmddhh_init2.longname = "Initial condition date/time in yyyymmddhh format"

    yyyymmddhh_fcst2 = rootgrp.createVariable('yyyymmddhh_fcst','i4',('time',))
    yyyymmddhh_fcst2.longname = "Forecast valid date/time in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    apcp_fcst = rootgrp.createVariable('apcp_fcst','f4',('time','ens','yf','xf',),
        zlib=True,least_significant_digit=2)  
    apcp_fcst.units = "mm" 
    apcp_fcst.long_name = "Ensemble precipitation forecast"
    apcp_fcst.valid_range = [0.,1000.]
    apcp_fcst.missing_value = np.array(-9999.99,dtype=np.float32)

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

