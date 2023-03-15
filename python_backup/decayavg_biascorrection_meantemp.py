
def decayavg_biascorrection_meantemp(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, date_list_forecast, alpha):
    
    """ apply modified decaying average bias correction 
        to forecasts.  Use the mean of all locations with temperatures 
        less than 1.5C different
    """
    
    import numpy as np
    from datetime import datetime
    import sys
    
    # -------------------------------------------------------------
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    # ---- sequentially loop through dates during the sample, 
    #      updating the previous day's bias correction 
    #      to the new days fcst vs. obs discrepancy.
    
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    ones = np.ones((nlats,nlons), dtype=np.int32)
    zeros = np.zeros((nlats,nlons), dtype=np.int32)
    for idate, date in enumerate(date_list_forecast):

        # ---- calculate a bias correction using the mean difference
        #      of F-O's for the points that have similar forecast temperatures.
        
        print ('processing date = ', date)
        forecast_2d = forecast_3d[idate,:,:] 
        analysis_2d = analyses_3d[idate,:,:]
        for i in range(nlons):
            for j in range(nlats):
                ftoday = forecast_2d[j,i]
                adiff = np.abs(forecast_2d - ftoday)
                a = np.where(adiff < 1.5, ones, zeros)
                fmean = np.sum(forecast_2d*a) / np.sum(ones)
                amean = np.sum(analysis_2d*a) / np.sum(ones)
                obsinc_2d[j,i] = fmean - amean
        #obsinc_2d[:,:] = forecast_3d[idate,:,:] - analyses_3d[idate,:,:]
            
        #if idate > 0: 
        #    obsinc_2d[:,:] = analyses_3d[idate,:,:] - \
        #        (forecast_3d[idate,:,:] - beta_3d[idate-1,:,:])
        #else:
        #    obsinc_2d[:,:] = analyses_3d[idate,:,:] - \
        #        forecast_3d[idate,:,:]
        
        if idate == 0:
            beta_3d[idate,:,:] = alpha*obsinc_2d[:,:]
        else:
            beta_3d[idate,:,:] = (1.-alpha)*beta_3d[idate-1,:,:] + alpha*obsinc_2d[:,:]

    return beta_3d