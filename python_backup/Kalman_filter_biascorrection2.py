
def Kalman_filter_biascorrection2(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, R, Bbeta,  \
    savefile, already_decay):
    
    """ apply decaying average Kalman filter bias correction 
        to forecasts.  Note the mix of some arrays shapes; 
        Kalman gain is shaped (nlats*nlons, nlats*nlons)
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    from reformat_gain_to_4d_f90 import reformat_gain_to_4d_f90
    
    # -------------------------------------------------------------
    
    def reformat_gain_to_4d(nlats, nlons, gain_2D):
    
        # ---- reform 2D Kalman gain into 4D-array.

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        gain_4D = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float64)
        ktr1 = 0
        for i1 in range(nlons):
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            for j1 in range(nlats):
                ktr2 = 0
                for i2 in range(nlons):
                    for j2 in range(nlats):
                        gain_4D[j1,i1,j2,i2] = gain_2D[ktr1,ktr2]
                        ktr2 = ktr2 + 1
                ktr1 = ktr1 + 1
        return gain_4D
    
    # -------------------------------------------------------------
    
    if already_decay == True:
        print ('   reading beta_3d from ', savefile) 
        inf = open(savefile, 'rb')
        beta_3d = cPickle.load(inf)
        inf.close()
        print ('   done reading')
    else:  
        # ---- define the Kalman gain per eq. (39) in Dick Dee
        #      Bias and Data Assimilation article, 
        #      https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.05.137
    
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('   max, min Bbeta = ',np.max(Bbeta), np.min(Bbeta))
        Bbeta_plus_R = R + Bbeta
        Bbeta_plus_R_inv = \
            np.linalg.inv(Bbeta_plus_R)
        Kalman_gain_beta = np.matmul(Bbeta, Bbeta_plus_R_inv)
        print ('   max, min Kalman_gain_beta  = ', \
            np.max(Kalman_gain_beta), np.min(Kalman_gain_beta))
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        Kalman_gain_beta_4d = reformat_gain_to_4d_f90(\
            Kalman_gain_beta, nlats, nlons)        
    
        # ---- sequentially loop through dates during the sample, 
        #      updating the previous day's bias correction 
        #      to the new days fcst vs. obs discrepancy.
    
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
                obsinc_2d[:,:] = analyses_3d[idate,:,:] - \
                    forecast_3d[idate,:,:]
        
            for i in range(nlons):
                for j in range(nlats):
                    ktr = nlats*j + i
                    
                    # ---- update the bias correction estimate, eq. 37 in Dee. 
                    
                    if idate > 0:
                        beta_3d[idate,j,i] = beta_3d[idate-1,j,i] - \
                        np.sum(Kalman_gain_beta_4d[j,i,:,:]*obsinc_2d[:,:]) 
        
                    else:
                        beta_3d[idate,j,i] = \
                            -np.sum(Kalman_gain_beta_4d[j,i,:,:]*\
                            obsinc_2d[:,:])
                                          
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('   ending decay avg Kalman_filter_biascorrection.py.'+\
            '  Current time = ', current_time)
        print ('   writing beta_3d to ', savefile) 
        ouf = open(savefile, 'wb')
        cPickle.dump(beta_3d, ouf)
        ouf.close()
        print ('   done writing')
    
    return beta_3d