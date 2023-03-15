def read_forecast_timeseries_GEFSv12_cloud(nlats, nlons, \
    cpath_forecast, date_list, clead):
    
    """
    read_forecast_timeseries_GEFSv12_cloud.py : read in the forecast data desired
    """
    
    import numpy as np
    import numpy.ma as ma
    import pygrib
    import os.path as path
    import sys
    import _pickle as cPickle

    ndates = len(date_list)
    for idate, date in enumerate(date_list):
 
        # ---- read the forecast information for bias corr.
         
        cyear = date[0:4]
        cmmdd = date[4:8]    
        infile = cpath_forecast + cyear + '/'+date+'_lead'+\
            clead+'_cldcover_conus_0.5deg_hour'+clead+'.cPick'
        fexist2 = path.exists(infile)
        if fexist2 == True:
            if cmmdd == '0101':
                forecast_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            inf = open(infile,'rb')
            forecast = cPickle.load(inf)
            nlats, nlons = np.shape(forecast)
            inf.close()
            forecast_3d[idate,:,:] = forecast[:,:]
        else:
            print ('Unable to read ', infile)
            forecast_3d[idate,:,:] = ma.masked
 
    return forecast_3d