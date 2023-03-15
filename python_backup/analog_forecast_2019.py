
def analog_forecast_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, clead, cpath_forecast, \
    cpath_era5):
    
    """ apply MOS forecast regression procedure in 2019. 
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    from dateutils import daterange, dateshift, dayofyear, splitdate
    
    
    cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
         
    # ---- read the climatology file.
    
    if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
        infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
    else:
        infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'
    inf = open(infile,'rb')
    climo_temps_estimated = cPickle.load(inf)
    inf.close()
         
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    fcollect = []
    acollect = []
    first = True
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    #for idate, date in enumerate(date_list_forecast[1:]):
    for idate, date in enumerate(date_list_forecast):

        cdd = date[6:8]
        cmm = date[4:6]
        cyearf = date[0:4]
        if first == True or cdd == '01':
            
            cmonth = cmonths[int(cmm)-1]
            first = False
            
            # ---- read the sorted forecasts and associated analyzed from file
            
            infile = cpath_forecast+'forecast_analyzed_sorted_'+cmonth+\
                '_lead='+clead+'.cPick'
            inf = open(infile,'rb')
            forecast_validdates_sorted = cPickle.load(inf)
            analysis_validdates_sorted = cPickle.load(inf)
            inf.close()
            nasamps, nalats, nalons = np.shape(forecast_validdates_sorted)
            
        # ---- find the closest forecast in the sorted data to today's forecast.
        #      form mean forecast and analyzed from nearest 25 samples.   Bias
        #      correction is then the mean(F) - mean(A)
        
        
        for jy in range(nlats):
            for ix in range(nlons):
        #for jy in range(nlats//2,nlats//2+1):
        #    for ix in range(nlons//2, nlons//2+1):
                if lsmask[jy,ix] == 1:
                    f = forecast_validdates_sorted[:,jy,ix]
                    a = analysis_validdates_sorted[:,jy,ix]
                    #print ('todays forecast = ', forecast_3d[idate,jy,ix])
                    #print ('f[0:-1:10] = ', f[0:-1:10])
                    #print ('a[0:-1:10] = ', a[0:-1:10])
                    idx = np.argmin(np.abs(f-forecast_3d[idate,jy,ix]))
                    idxmin = np.max([idx-21,0])
                    idxmax = np.min([nasamps,idx+21])
                    #print ('idx, idxmin, idxmax = ', idx, idxmin, idxmax)
                    #print ('f[idxmin:idxmax] = ', f[idxmin:idxmax])
                    #print ('a[idxmin:idxmax] = ', a[idxmin:idxmax])
                    #print ('fmean, amean = ', np.mean(f[idxmin:idxmax]),  np.mean(a[idxmin:idxmax]))
                    beta_3d[idate,jy,ix] = np.mean(f[idxmin:idxmax]) - \
                        np.mean(a[idxmin:idxmax])
                    #print ('beta_3d = ', beta_3d[idate,jy,ix])
                    #sys.exit()
                
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    return beta_3d