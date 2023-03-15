
def lowess_forecast_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, clead, cpath_forecast, \
    cpath_era5):
    
    """ apply lowess regression procedure in 2019. 
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    from dateutils import daterange, dateshift, dayofyear, splitdate
    from implement_lowess import implement_lowess
    
    cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
     
    #print ('begin lowess forecast 2019')    
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
        print (idate, date)
        cdd = date[6:8]
        cmm = date[4:6]
        cyearf = date[0:4]
        if first == True or cdd == '01':
            #print (first, cdd)
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
            
            # ----- lowess regression.   return an ordered set of regular forecast
            #       values and their predicted analyzed values
            
            fcstvals = np.zeros((26, nlats, nlons), dtype=np.float64)
            predicted_anal_vals = np.zeros((26, nlats, nlons), dtype=np.float64)
            
            #print ('forecast_validdates_sorted = ',forecast_validdates_sorted[:,1,0])
            #print ('analysis_validdates_sorted = ',analysis_validdates_sorted[:,1,0])
            fcstvals, predicted_anal_vals = implement_lowess(\
                forecast_validdates_sorted, analysis_validdates_sorted, \
                lsmask, nasamps, nalats, nalons, fcstvals, predicted_anal_vals)
                
        # ---- lowess regression returns predicted analyzed at regular set of 
        #      forecast values.   Use the difference between forecast and
        #      predicted as the bias correction
        
        for jy in range(nlats):
            for ix in range(nlons):
        #for jy in range(nlats//2, nlats//2+1):
        #    for ix in range(nlons//2, nlons//2+1):
                if lsmask[jy,ix] == 1:
                    f = fcstvals[:,jy,ix]
                    a = predicted_anal_vals[:,jy,ix]
                    idx = np.argmin(np.abs(f-forecast_3d[idate,jy,ix]))
                    beta_3d[idate,jy,ix] = f[idx] - a[idx]
                    
                    #if jy == 1 and ix == 0:
                        #print ('lsmask[1,0] = ', lsmask[1,0])
                        #print ('f = ', f)
                        #print ('a = ', a)
                        #print ('forecast_3d[today] = ', forecast_3d[idate,jy,ix])
                        #print ('beta, f, a = ', beta_3d[idate,jy,ix], f[idx], a[idx])
                        #if f[idx] - a[idx] == 0.0:
                        #    print ('processing grid point jy, ix = ', jy, ix)
                        #    print ('f = ', f)
                        #    print ('a = ', a)
                        #    sys.exit()

                    #sys.exit()
                    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    return beta_3d