
def Kalman_filter_biascorrection_2018gain(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, Kalman_gain_beta_4D):
    
    """ apply Kalman filter bias correction to forecasts.  
    """
    import numpy as np
    from datetime import datetime
    import sys
         
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    beta_2d = np.zeros((nlats, nlons), dtype=np.float64)
    for idate, date in enumerate(date_list_forecast[1:]):

        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)
            
        if idate > 0: 
            obsinc_2d[:,:] = analyses_3d[idate,:,:] - \
                (forecast_3d[idate,:,:] - beta_3d[idate-1,:,:])
        else:
            obsinc_2d[:,:] = analyses_3d[idate,:,:] - forecast_3d[idate,:,:]
        
        for i in range(nlons):
            for j in range(nlats):
                ktr = nlats*j + i
                    
                # ---- update the bias correction estimate, eq. 37 in Dee. 
                    
                if idate > 0:
                    beta_3d[idate,j,i] = beta_3d[idate-1,j,i] - \
                        np.sum(Kalman_gain_beta_4D[j,i,:,:]*obsinc_2d[:,:]) 
                    e = beta_3d[idate-1,j,i]
                else:
                    beta_3d[idate,j,i] = -np.sum(Kalman_gain_beta_4D[j,i,:,:]*obsinc_2d[:,:])
                    e = 0.0
                                          
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    #print ('ending Kalman_filter_biascorrection.py.  Current time = ', current_time)
    
    return beta_3d