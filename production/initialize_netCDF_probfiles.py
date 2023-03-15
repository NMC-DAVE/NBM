def initialize_netCDF_probfiles(outfile, ncname, ny, nx, \
    lats_in, lons_in, nthresholds, thresholds):
        
    """ initialize the netCDF files for writing probability output"""

    from netCDF4 import Dataset
    import sys
    import numpy as np


    print ('initializing netCDF file ',outfile)
    if ncname == 'nc_fullfield':
        nc_fullfield = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
        ncout = nc_fullfield
    elif ncname == 'nc_thinned':
        nc_thinned = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
        ncout = nc_thinned
    else:
        print ('invalid ncname in initalize_netCDF_probfiles = ', ncname)
        sys.exit()
        
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
    
    threshold = ncout.createDimension('threshold',nthresholds)
    thresholdv = ncout.createVariable('thresholdv','i4',('threshold',))
    thresholdv.units = "mm"

    yyyymmddhh_init = ncout.createVariable('yyyymmddhh_init','i4',('sample',))
    yyyymmddhh_init.longname = "Initial condition date/time in yyyymmddhh format"

    yyyymmddhh_fcst = ncout.createVariable('yyyymmddhh_fcst','i4',('sample',))
    yyyymmddhh_fcst.longname = "Forecast valid date/time in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid
    
    probability_raw = ncout.createVariable('probability_raw','f4',\
        ('sample','threshold','yf','xf',),zlib=True,least_significant_digit=3)
    probability_raw.units = "n/a"
    probability_raw.long_name = 'Probability of precipitation '+\
        'exceeding threshold amount in raw ens.'
    probability_raw.valid_range = [0.,1.]
    probability_raw.missing_value = np.array(-99.99,dtype=np.float32)
    
    probability_qmapped = ncout.createVariable('probability_qmapped',\
        'f4',('sample','threshold','yf','xf',),
        zlib=True,least_significant_digit=3)
    probability_qmapped.units = "n/a"
    probability_qmapped.long_name = 'Probability of precipitation '+\
        'exceeding threshold amount in qmapped ensemble'
    probability_qmapped.valid_range = [0.,1.]
    probability_qmapped.missing_value = np.array(-99.99,dtype=np.float32)
    
    probability_qmapped_weighted = ncout.createVariable('probability_qmapped_weighted',\
        'f4',('sample','threshold','yf','xf',),
        zlib=True,least_significant_digit=3)
    probability_qmapped_weighted.units = "n/a"
    probability_qmapped_weighted.long_name = 'Probability of quantile mapped and '+\
        'closest-member histogram weighted precipitation '+\
        'exceeding threshold amount in qmapped ensemble'
    probability_qmapped_weighted.valid_range = [0.,1.]
    probability_qmapped_weighted.missing_value = np.array(-99.99,dtype=np.float32)
    
    probability_qmapped_weighted_dressed = \
        ncout.createVariable('probability_qmapped_weighted_dressed',\
        'f4',('sample','threshold','yf','xf',),
        zlib=True,least_significant_digit=3)
    probability_qmapped_weighted_dressed.units = "n/a"
    probability_qmapped_weighted_dressed.long_name = 'Probability of quantile mapped and '+\
        'closest-member histogram weighted and Gaussian-dressed precipitation '+\
        'exceeding threshold amount in qmapped ensemble'
    probability_qmapped_weighted_dressed.valid_range = [0.,1.]
    probability_qmapped_weighted_dressed.missing_value = np.array(-99.99,dtype=np.float32)
    
    ensmean_raw = ncout.createVariable('ensmean_raw','f4',\
        ('sample','yf','xf',),zlib=True,least_significant_digit=3)
    ensmean_raw.units = "mm"
    ensmean_raw.long_name = 'Raw ensemble-mean precicipitation amount (mm)'
    ensmean_raw.valid_range = [0.,200.]
    ensmean_raw.missing_value = np.array(-99.99,dtype=np.float32)
    
    ensmean_qmapped = ncout.createVariable('ensmean_qmapped',\
        'f4',('sample','yf','xf',),
        zlib=True,least_significant_digit=3)
    ensmean_qmapped.units = "mm"
    ensmean_qmapped.long_name = 'Quantile-mapped ensemble-mean '+\
        'precicipitation amount (mm)'
    ensmean_qmapped.valid_range = [0.,200.]
    ensmean_qmapped.missing_value = np.array(-99.99,dtype=np.float32)
    
    # ---- metadata

    ncout.title = 'NDFD CONUS gridded probabilities of precipitation exceeding '+\
        'various thresholds, weighted and quantile mapped'
    ncout.history = "GEFSv12 implemented at NCEP/EMC Sep 2020"
    ncout.institution =  "NCEP/EMC and PSL"
    ncout.platform = ""
    ncout.references = ""

    # ---- initialize

    xvf[:] = np.arange(nx)
    yvf[:] = np.arange(ny)
    lons_out[:] = lons_in[:,:]
    lats_out[:] = lats_in[:,:]
    thresholdv[:] = thresholds[:]
    
    if ncname == 'nc_fullfield':
        nc_fullfield.close()
    elif ncname == 'nc_thinned':
        nc_thinned.close() 
    
    istat = 0
    
    return istat
    