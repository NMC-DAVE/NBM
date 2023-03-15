
def Kalman_filter_biascorrection_savegain2(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, R, Bbeta, b_beta_var, gainfile, already):
    
    """ apply Kalman filter-like bias correction to forecasts.  Note
        the mix of some arrays shapes; Kalman gain is shaped
        (nlats*nlons, nlats*nlons)
    """
    
    import numpy as np
    from datetime import datetime
    import sys
    from reformat_gain_to_4d_f90 import reformat_gain_to_4d_f90
    from update_beta import update_beta
    import _pickle as cPickle
    
    # -------------------------------------------------------------
    
    def reformat_gain_to_4d(nlats, nlons, gain_2D):
    
        # ---- reform 2D Kalman gain into 4D-array.

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        #print ('Reforming into 4D Kalman gain matrix. Current time = ', current_time)
        gain_4D = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float64)
        ktr1 = 0
        for i1 in range(nlons):
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            #print ('processing i1 = ',i1,' of ',nlons,'. Current time = ',current_time)
            for j1 in range(nlats):
                ktr2 = 0
                for i2 in range(nlons):
                    for j2 in range(nlats):
                        gain_4D[j1,i1,j2,i2] = gain_2D[ktr1,ktr2]
                        ktr2 = ktr2 + 1
                ktr1 = ktr1 + 1
        return gain_4D
    
    # -------------------------------------------------------------
    
    #print (update_beta.__doc__)
    if already == False:   
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        Bbeta_plus_R = R + Bbeta*b_beta_var*0.01
        Bbeta_plus_R_inv = np.linalg.inv(Bbeta_plus_R)
        Kalman_gain_beta = np.matmul(Bbeta*b_beta_var*0.01, Bbeta_plus_R_inv)
        print ('   max, min Bbeta = ',np.max(Bbeta), np.min(Bbeta))
        print ('   max, min Kalman_gain_beta  = ', \
            np.max(Kalman_gain_beta), np.min(Kalman_gain_beta))
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        Kalman_gain_beta_4d = reformat_gain_to_4d_f90(Kalman_gain_beta, nlats, nlons)  
    
        print ('   writing Kalman gain to ', gainfile) 
        ouf = open(gainfile, 'wb')
        cPickle.dump(Kalman_gain_beta_4d, ouf)
        ouf.close()
        print ('   done writing')
    else:
        print ('   reading Kalman gain from ', gainfile) 
        inf = open(gainfile, 'rb')
        Kalman_gain_beta_4d = cPickle.load(inf)
        inf.close()
        print ('   done reading')
         
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    beta_2d = np.zeros((nlats, nlons), dtype=np.float64)
    for idate, date in enumerate(date_list_forecast[1:]):

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
                    
                # ---- update the bias correction estimate 
                    
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
    #print ('ending Kalman_filter_biascorrection.py.  Current time = ', current_time)
    
    return beta_3d