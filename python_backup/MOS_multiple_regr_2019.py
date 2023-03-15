
def MOS_multiple_regr_2019(npts, nlats, nlons, \
    analyses_3d, forecast_3d, beta_decay_3d, beta_qmap_3d, beta_3d, \
    lsmask, date_list_forecast, clead, cpath_forecast, \
    cpath_era5):
    
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
         
    # ---- read the climatology file.
    
    if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
        infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
    else:
        infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'
    inf = open(infile,'rb')
    climo_temps_estimated = cPickle.load(inf)
    inf.close()
         
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
            
            infile = cpath_forecast+'MOS_multiple_regr_'+cmonth+\
                '_lead='+clead+'.cPick'
            print (idate, date, infile)
            inf = open(infile,'rb')
            regr_coefs = cPickle.load(inf)
            inf.close()
            intercept = regr_coefs[0,:,:]
            slope = regr_coefs[1,:,:]
            beta_slope = regr_coefs[2,:,:]
            qmap_slope = regr_coefs[3,:,:]
            #qm_beta_interact = regr_coefs[4,:,:]
            #f_beta_interact = regr_coefs[5,:,:]
            #print ('date, int sl beta qmap = ', \
            #    idate, date, intercept[nlats//2, nlons//2], \
            #    slope[nlats//2, nlons//2], \
            #    beta_slope[nlats//2, nlons//2],\
            #    qmap_slope[nlats//2, nlons//2])
            
        # ---- determine the julian day

        imm = int(cmm)
        idd = int(cdd)
        iyear_full = int(cyearf)
        julday = dayofyear(iyear_full, imm, idd) - 1
        if julday > 364: julday = 364 
         
        #ftoday = forecast_3d[idate,:,:] - climo_temps_estimated[julday,:,:]
        #atoday = analyses_3d[idate,:,:] - climo_temps_estimated[julday,:,:]
        ftoday = forecast_3d[idate,:,:] 
        atoday = analyses_3d[idate,:,:] 
        
        if idate - iskip >= 0:
            beta_today = beta_decay_3d[idate-iskip,:,:]
        else:
            beta_today[:,:] = 0.0
        qmap_today = beta_qmap_3d[idate,:,:]

        idxmax = np.argmax(qmap_today*qmap_slope)
        idxmin = np.argmin(qmap_today*qmap_slope)
        regressed = intercept[:,:] + slope[:,:]*ftoday[:,:] + \
            beta_slope[:,:]*beta_today[:,:] + \
            qmap_slope[:,:]*qmap_today[:,:] 
        beta_3d[idate,:,:] = lsmask[:,:]*(ftoday[:,:]-regressed[:,:])
        
        
        
        #print (idate, ftoday[nlats//2, nlons//2], qmap_today[nlats//2, nlons//2], \
        #    beta_today[nlats//2, nlons//2], regressed[nlats//2, nlons//2])
            
        #diff = analyses_3d[idate,59,72] - \
        #    (forecast_3d[idate,59,72] - beta_3d[idate,59,72])
        #if np.abs(diff) > 10.0:
        #    print (idate,date,analyses_3d[idate,59,72], forecast_3d[idate,59,72], beta_3d[idate,59,72] )
        #    print ('idate, fprime, f, cl = ', idate, ftoday[59,72], forecast_3d[idate,59,72], \
        #        climo_temps_estimated[julday,59,72])
        #    print ('    beta_decay, qmap_today = ', beta_today[59,72], qmap_today[59,72])
        #    print ('    i,s, bs, qs = ', intercept[59,72], slope[59,72], beta_slope[59,72], qmap_slope[59,72])
        
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    
 
        
    return beta_3d