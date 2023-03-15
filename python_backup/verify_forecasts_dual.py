def verify_forecasts_dual(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, bias_corr_3d, bias_qmap, lsmask, \
    iskip, ftype, statsfile):
    
    """ 
    verify_forecasts_dual: over the valid samples, evaluate
    the bias-corrected forecasts
    """
    
    import numpy as np
    import numpy.ma as ma
    import _pickle as cPickle
    
    
    rmse = np.float64(0.0)
    mae = np.float64(0.0)
    bia = np.float64(0.0)
    rmse_byday = ma.zeros((ndates), dtype=np.float64)
    mae_byday = ma.zeros((ndates), dtype=np.float64)
    bia_byday = ma.zeros((ndates), dtype=np.float64)
    difference = ma.zeros((nlats, nlons), dtype=np.float64)
    for idate in range(iskip, ndates-iskip):
            
        bias = (bias_corr_3d[idate-iskip,:,:] + bias_qmap[idate,:,:])/2.
        difference[:,:] = analyses_3d[idate,:,:] - \
            (forecast_3d[idate,:,:] - bias[:,:])
        rmse = rmse + np.sum(lsmask*difference**2)
        mae = mae + np.sum(np.abs(lsmask*difference))
        bia = bia + np.sum(lsmask*difference)
        
        rmse_byday[idate] = np.sqrt(np.sum(lsmask*difference**2)) / float(np.sum(lsmask))
        mae_byday[idate] = np.sum(np.abs(lsmask*difference)) / float(np.sum(lsmask))
        bia_byday[idate] = np.sum(lsmask*difference) / float(np.sum(lsmask))
        
    nsamps = (ndates-(2*iskip))*np.sum(lsmask)
    bia = bia / float(nsamps)
    mae = mae / float(nsamps)
    rmse = np.sqrt(rmse/nsamps)
    
    # --- write to text file
    
    ouf = open(statsfile, 'w')
    print (bia, mae, rmse, file=ouf)
    ouf.close()
    
    statsfile2 = statsfile[0:-4]+'.cPick'
    print ('writing to ', statsfile2)
    ouf = open(statsfile2, 'wb')
    cPickle.dump(rmse_byday, ouf)
    cPickle.dump(bia_byday, ouf)
    cPickle.dump(mae_byday, ouf)
    ouf.close()
    
    return rmse, bia, mae 