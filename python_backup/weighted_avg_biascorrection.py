
def weighted_avg_biascorrection(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, Bbeta, alpha):
    
    """ apply weighted average bias correction, a spatial 
        extension of decaying average
    """
    import numpy as np
    from datetime import datetime
    import sys
    
    # -------------------------------------------------------------
    
    # ---- sequentially loop through dates during the sample, 
    #      updating the previous day's bias correction 
    #      to the new days fcst vs. obs discrepancy.
    
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    beta_2d = np.zeros((nlats, nlons), dtype=np.float64)
    for idate, date in enumerate(date_list_forecast):

        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)
            
        obsinc_2d[:,:] = forecast_3d[idate,:,:] - analyses_3d[idate,:,:] 
        
        for i in range(nlons):
            for j in range(nlats):
                ktr = nlats*j + i
                    
                # ---- weight the observation increment by correlation function
                #      for this point.   Then adjust bias with b_beta_var kludge
                #      factor * weighted obs increment
                    
                numer = np.sum(Bbeta[j,i,:,:]*obsinc_2d[:,:])
                denom = np.sum(Bbeta[j,i,:,:])
                if idate > 0:
                    beta_3d[idate,j,i] = (1.-alpha)*beta_3d[idate-1,j,i] + \
                        alpha*numer/denom 
                else:
                    beta_3d[idate,j,i] = numer/denom
                                          
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('   ending weighted_avg_biascorrection.py.'+\
        '  Current time = ', current_time)
    
    return beta_3d