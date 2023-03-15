def read_forecast_timeseries(cpath_forecast, date_list, clead, cvariable):
    
    """
    read_forecast_timeseries.py : read in the ERA-5 reanalysis data desired
    return 3D array with data as well as associated lat/lon
    """
    
    import numpy as np
    import numpy.ma as ma
    import pygrib
    import os.path as path
    import sys

    ndates = len(date_list)
    read_latlons = False
    for idate, date in enumerate(date_list):
 
        # ---- read the control forecast at this lead time and initial date

        infile = cpath_forecast + cvariable+'_'+date+'_f'+clead+'.grib2'  
        fexist = path.exists(infile)
        #print (infile, fexist)
        if fexist == True:
            grbfile = pygrib.open(infile) 
            grb = grbfile.select()[0] 
            forecast = grb.values
            if idate == 0:
                nlats, nlons = np.shape(forecast)
                forecast_3d = np.zeros((ndates,nlats,nlons))
                lats, lons = grb.latlons()
                lats = np.flipud(lats)
                #print ('forecast lats[:,0] = ', lats[:,0])
            grbfile.close()
            forecast_3d[idate,:,:] = np.flipud(forecast[:,:])
        else:
            print ('Unable to read ', infile)
            forecast_3d[idate,:,:] = ma.masked 
 
    return forecast_3d, lats, lons