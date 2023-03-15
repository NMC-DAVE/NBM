
def Kalman_filter_biascorr_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, clead, cpath_gain):
    
    """ apply Kalman filter bias correction to forecasts in 2019. 
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    
    # -------------------------------------------------------------
    
    # ---- based roughly on what had the lowest error in 2018, set the Kalman
    #      localization radii used
    
    if clead == '24':
        flocal_warm = '400.0'
        blocal_warm = '600.0' # '1200.0'
        flocal_cold = '400.0'
        blocal_cold = '600.0' # '1200.0'
    elif clead == '48':
        flocal_warm = '400.0'
        blocal_warm = '600.0' # ''1200.0'
        flocal_cold = '400.0'
        blocal_cold = '600.0' # '1200.0'
    elif clead == '72':
        flocal_warm = '400.0'
        blocal_warm = '1200.0'
        flocal_cold = '600.0'
        blocal_cold = '1200.0'
    elif clead == '96':
        flocal_warm = '400.0'
        blocal_warm = '1200.0'
        flocal_cold = '600.0'
        blocal_cold = '1200.0'
    elif clead == '120':
        flocal_warm = '200.0'
        blocal_warm = '1200.0'
        flocal_cold = '600.0'
        blocal_cold = '1200.0'   
    
    # ---- define the Kalman gain per eq. (39) in Dick Dee
    #      Bias and Data Assimilation article, 
    #      https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.05.137
   
    gainfile_cold = cpath_gain + '2018_KFgain_flocal'+\
        flocal_cold+'_blocal'+blocal_cold+\
        '_2018_cold_lead'+clead+'.cPick'
    gainfile_warm = cpath_gain + '2018_KFgain_flocal'+\
        flocal_warm+'_blocal'+blocal_warm+\
        '_2018_warm_lead'+clead+'.cPick'

    print ('reading cold_season Kalman gain from ', gainfile_cold) 
    inf = open(gainfile_cold, 'rb')
    Kalman_gain_beta_4d_cold = cPickle.load(inf)
    inf.close()
    print ('done reading')
    
    print ('reading warm_season Kalman gain from ', gainfile_warm) 
    inf = open(gainfile_warm, 'rb')
    Kalman_gain_beta_4d_warm = cPickle.load(inf)
    inf.close()
    print ('done reading')
         
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    beta_2d = np.zeros((nlats, nlons), dtype=np.float64)
    for idate, date in enumerate(date_list_forecast[1:]):

        # ---- determine which gain file to use
        
        if idate == 0 : #'2019010100':
            Kalman_gain_beta_4d = Kalman_gain_beta_4d_cold
        elif date == '2019040100':
            Kalman_gain_beta_4d = Kalman_gain_beta_4d_warm
        elif date == '2019100100':
            Kalman_gain_beta_4d = Kalman_gain_beta_4d_cold

        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        #print ('processing date, current time = ', date, current_time ) 
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
                        np.sum(Kalman_gain_beta_4d[j,i,:,:]*obsinc_2d[:,:]) 
                else:
                    beta_3d[idate,j,i] = -np.sum(Kalman_gain_beta_4d[j,i,:,:]* \
                        obsinc_2d[:,:])
                
        print ('idate, date, max, min beta_3d = ', idate, date, \
            np.max(beta_3d[idate,:,:]), np.min(beta_3d[idate,:,:]) )
                                          
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('finishing Kalman_filter_biascorr_2019. ', current_time)
    
    return beta_3d