
def KFgain_GEFSv12_biascorr_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, clead, cpath_gain):
    
    """ apply Kalman filter bias correction to forecasts in 2019. 
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    
    cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
         
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    first = True
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    for idate, date in enumerate(date_list_forecast[1:]):

        cdd = date[6:8]
        cmm = date[4:6]
        if first == True or cdd == '01':
            
            cmonth = cmonths[int(cmm)-1]
            first = False
            
            # ---- define the Kalman gain per eq. (39) in Dick Dee
            #      Bias and Data Assimilation article, 
            #      https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.05.137
    
            gain_infile = cpath_gain + 'GEFSv12_KFgain_'+cmonth+'_lead'+clead+'.cPick'
            inf = open(gain_infile, 'rb')
            Kalman_gain_beta_4d = cPickle.load(inf)
            inf.close()

        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        if idate > 0: 
            obsinc_2d[:,:] = analyses_3d[idate,:,:] - \
                (forecast_3d[idate,:,:] - beta_3d[idate-1,:,:])
        else:
            obsinc_2d[:,:] = analyses_3d[idate,:,:] - forecast_3d[idate,:,:]
        
        for i in range(nlons):
            for j in range(nlats):
                    
                # ---- update the bias correction estimate, eq. 37 in Dee. 
                    
                #if lsmask[j,i] == 1:
                #    if idate > 0:
                #        beta_3d[idate,j,i] = beta_3d[idate-1,j,i] - \
                #            np.sum(Kalman_gain_beta_4d[j,i,:,:]*obsinc_2d[:,:]) 
                #   else:
                #        beta_3d[idate,j,i] = -np.sum(Kalman_gain_beta_4d[j,i,:,:]* \
                #            obsinc_2d[:,:])
                #else:
                #    beta_3d[idate,j,i] = 0.0
                    
                    
                if idate > 0:
                    beta_3d[idate,j,i] = beta_3d[idate-1,j,i] - \
                        np.sum(Kalman_gain_beta_4d[j,i,:,:]*obsinc_2d[:,:]) 
                else:
                    beta_3d[idate,j,i] = -np.sum(Kalman_gain_beta_4d[j,i,:,:]* \
                        obsinc_2d[:,:])
                                          
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('finishing KFgain_GEFSv12_biascorr_2019. ', current_time)
    
    return beta_3d