def verify_forecasts_2000_2018(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, bias_corr_3d, lsmask, \
    lons, lats, iskip, statsfile):
    
    """ 
    verify_forecasts_2000_2018: over the valid samples, evaluate
    the bias-corrected forecasts
    """
    
    import numpy as np
    import numpy.ma as ma
    import _pickle as cPickle
    import sys
    
    ibracket = int(clead) // 24
    rmse = np.float64(0.0)
    mae = np.float64(0.0)
    bia = np.float64(0.0)
    
    difference = ma.zeros((nlats, nlons), dtype=np.float64)
    
    nactdates = 0
    for idate in range(ibracket, ndates-ibracket):
           
        nactdates = nactdates+1 
        difference[:,:] = analyses_3d[idate,:,:] - \
            (forecast_3d[idate,:,:] - bias_corr_3d[idate-iskip,:,:])
        
        #print ('idate, rmse = ', idate, np.sqrt(np.sum(lsmask*difference**2)/np.sum(lsmask)))
        rmse = rmse + np.sum(lsmask*difference**2)
        mae = mae + np.sum(np.abs(lsmask*difference))
        bia = bia + np.sum(lsmask*difference) 
        
    nsamps = nactdates*np.sum(lsmask)
    bia = bia / float(nsamps)
    mae = mae / float(nsamps)
    rmse = np.sqrt(rmse/nsamps)
    
    # --- write to text file
    
    ouf = open(statsfile, 'w')
    print (bia, mae, rmse, file=ouf)
    ouf.close()
    
    return rmse, bia, mae