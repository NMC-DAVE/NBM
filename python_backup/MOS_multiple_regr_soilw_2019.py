def MOS_multiple_regr_soilw_2019(npts, nlats, nlons, \
    analyses_3d, forecast_3d, forecast_3d_soilw, forecast_3d_cloud, \
    beta_3d, lsmask, date_list_forecast, clead, cpath_forecast):
    
    """ apply MOS multiple regression procedure in 2019. 
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    from dateutils import daterange, dateshift, dayofyear, splitdate
    
    cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    iskip = int(clead) // 24
    if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
        iskip = iskip+1
         

    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    first = True
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    beta_today = np.zeros((nlats, nlons), dtype=np.float64)


    #for idate, date in enumerate(date_list_forecast[1:]):
    for idate, date in enumerate(date_list_forecast):

        cdd = date[6:8]
        cmm = date[4:6]
        cyearf = date[0:4]
        if first == True or cdd == '01':
            
            cmonth = cmonths[int(cmm)-1]
            first = False
            
            # ---- read the MOS regression coefficient data
            
            infile = cpath_forecast+'MOS_soil_cloud_'+cmonth+\
                '_lead='+clead+'.cPick'
            print (idate, date, infile)
            inf = open(infile,'rb')
            regr_coefs = cPickle.load(inf)
            inf.close()
            intercept = regr_coefs[0,:,:]
            slope = regr_coefs[1,:,:]
            soilw_slope = regr_coefs[2,:,:]
            cloud_slope = regr_coefs[3,:,:]
            soil_interact = regr_coefs[4,:,:]
            cloud_interact = regr_coefs[5,:,:]   
            cloud_soil_interact = regr_coefs[6,:,:]  
            
        # ---- determine the julian day

        imm = int(cmm)
        idd = int(cdd)
        iyear_full = int(cyearf)
        julday = dayofyear(iyear_full, imm, idd) - 1
        if julday > 364: julday = 364 
         
        ftoday = forecast_3d[idate,:,:] 
        ftoday_soilw = forecast_3d_soilw[idate,:,:] 
        ftoday_cloud = forecast_3d_cloud[idate,:,:] 
        atoday = analyses_3d[idate,:,:] 
    
        regressed = intercept[:,:] + slope[:,:]*ftoday[:,:] + \
            soilw_slope[:,:]*ftoday_soilw[:,:] + \
            cloud_slope[:,:]*ftoday_cloud[:,:] + \
            soil_interact[:,:]*ftoday[:,:]*ftoday_soilw[:,:] + \
            cloud_interact[:,:]*ftoday[:,:]*ftoday_cloud[:,:] + \
            cloud_soil_interact[:,:]*ftoday_soilw[:,:]*ftoday_cloud[:,:] 
        beta_3d[idate,:,:] = lsmask[:,:]*(ftoday[:,:]-regressed[:,:])
        
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    return beta_3d