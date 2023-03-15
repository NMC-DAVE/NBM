
def simple_decayavg_biascorrection(nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, alpha):
    
    """ apply decaying average bias correction
    """
    import numpy as np
    from datetime import datetime
    import sys
        
    # -------------------------------------------------------------
    
    
    # ---- define the Kalman gain per eq. (39) in Dick Dee
    #      Bias and Data Assimilation article, 
    #      https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.05.137
        

    
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    ndates = int(len(date_list_forecast))
    for idate, date in enumerate(date_list_forecast[1:]):

        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)
            
        if idate > 0: 
            beta_3d[idate,:,:] = (1-alpha)*beta_3d[idate-1,:,:] + \
                alpha*(forecast_3d[idate,:,:] - analyses_3d[idate,:,:])
    
    return beta_3d