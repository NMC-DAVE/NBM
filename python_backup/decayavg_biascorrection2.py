
def decayavg_biascorrection2(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, alpha):
    
    """ apply decaying average bias correction 
        to forecasts.  
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
    beta_2d = np.zeros((nlats, nlons), dtype=np.float64)
    for idate, date in enumerate(date_list_forecast):

        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)
        
        obsinc_2d[:,:] = lsmask[:,:]*(forecast_3d[idate,:,:] - analyses_3d[idate,:,:])
        
        if idate == 0:
            beta_3d[idate,:,:] = alpha*obsinc_2d[:,:]
        else:
            beta_3d[idate,:,:] = (1.-alpha)*beta_3d[idate-1,:,:] + \
                alpha*obsinc_2d[:,:]

    return beta_3d