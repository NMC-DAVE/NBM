def verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, bias_corr_3d, lsmask, \
    lons, lats, climo_temps_estimated, iskip, ftype, statsfile):
    
    """ 
    verify_forecasts: over the valid samples, evaluate
    the bias-corrected forecasts
    """
    
    import numpy as np
    import numpy.ma as ma
    import _pickle as cPickle
    import sys
    
    ibracket = int(clead) // 24
    
    if ftype == 'quantile_mapping_errorstats_2019_' or ftype == 'MOS_GEFSv12_2019_' \
    or ftype == 'MOS_mvr_GEFSv12_2019_':
        iskip = 0
    else:
        iskip = int(clead) // 24 
        if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
            iskip = iskip+1
        
    #print ('iskip = ', iskip,' ftype = ', ftype)
    rmse = np.float64(0.0)
    mae = np.float64(0.0)
    bia = np.float64(0.0)
    rmse_byday = ma.zeros((ndates), dtype=np.float64)
    mae_byday = ma.zeros((ndates), dtype=np.float64)
    bia_byday = ma.zeros((ndates), dtype=np.float64)
    
    rmse_bypoint = ma.zeros((nlats,nlons), dtype=np.float64)
    mae_bypoint = ma.zeros((nlats,nlons), dtype=np.float64)
    bia_bypoint = ma.zeros((nlats,nlons), dtype=np.float64)
    
    difference = ma.zeros((nlats, nlons), dtype=np.float64)
    forecast_3d_bcorr_prime = ma.zeros((ndates, nlats, nlons), dtype=np.float64)
    analysis_3d_prime = ma.zeros((ndates, nlats, nlons), dtype=np.float64)
    
    corr_winter = ma.zeros((nlats, nlons), dtype=np.float64)
    corr_spring = ma.zeros((nlats, nlons), dtype=np.float64)
    corr_summer = ma.zeros((nlats, nlons), dtype=np.float64)
    corr_fall   = ma.zeros((nlats, nlons), dtype=np.float64)
    std_normalized_winter = ma.zeros((nlats, nlons), dtype=np.float64)
    std_normalized_spring = ma.zeros((nlats, nlons), dtype=np.float64)
    std_normalized_summer = ma.zeros((nlats, nlons), dtype=np.float64)
    std_normalized_fall   = ma.zeros((nlats, nlons), dtype=np.float64)
    std_obs_winter = ma.zeros((nlats, nlons), dtype=np.float64)
    std_obs_spring = ma.zeros((nlats, nlons), dtype=np.float64)
    std_obs_summer = ma.zeros((nlats, nlons), dtype=np.float64)
    std_obs_fall   = ma.zeros((nlats, nlons), dtype=np.float64)
    
    forecast_3d_bcorr_prime = ma.zeros((ndates, nlats, nlons), dtype=np.float64)
    analysis_3d_prime = ma.zeros((ndates, nlats, nlons), dtype=np.float64)
    
    forecast_3d_bcorr_prime.mask = True
    analysis_3d_prime.mask = True
    
    nactdates = 0
    for idate in range(ibracket, ndates):
           
        nactdates = nactdates+1 
        
        difference[:,:] = analyses_3d[idate,:,:] - \
            (forecast_3d[idate,:,:] - bias_corr_3d[idate-iskip,:,:])
        #forecast_3d_bcorr_prime[idate,:,:] = (forecast_3d[idate,:,:] - \
        #    bias_corr_3d[idate-iskip,:,:]) - climo_temps_estimated[idate,:,:]
        #analysis_3d_prime[idate,:,:] = analyses_3d[idate,:,:] - \
        #    climo_temps_estimated[idate,:,:]
        forecast_3d_bcorr_prime[idate,:,:] = (forecast_3d[idate,:,:] - \
            bias_corr_3d[idate-iskip,:,:])
        analysis_3d_prime[idate,:,:] = analyses_3d[idate,:,:] 
        
        #print ('  idate, rmse = ', idate, \
        #    np.sqrt(np.sum(lsmask*difference**2)/np.sum(lsmask)))
        rmse = rmse + np.sum(lsmask*difference**2)
        mae = mae + np.sum(np.abs(lsmask*difference))
        bia = bia + np.sum(lsmask*difference)
        
        rmse_byday[idate] = \
            np.sqrt(np.sum(lsmask*difference**2) / float(np.sum(lsmask)))
        mae_byday[idate] = \
            np.sum(np.abs(lsmask*difference))  / float(np.sum(lsmask))
        bia_byday[idate] = \
            np.sum(lsmask*difference) / float(np.sum(lsmask))
            
        rmse_bypoint[:,:] = rmse_bypoint[:,:] + \
            lsmask*difference**2 
        mae_bypoint[:,:] = mae_bypoint[:,:] + \
            np.sum(np.abs(lsmask*difference))  
        bia_bypoint[:,:] = bia_bypoint[:,:]  +\
            np.sum(lsmask*difference)   
        
        #if ftype == 'MOS' or ftype == 'MOS_multiple_regression':
        #    print ('idate, difference = ',idate, difference[59,72])    
        
    nsamps = nactdates*np.sum(lsmask)
    bia = bia / float(nsamps)
    mae = mae / float(nsamps)
    rmse = np.sqrt(rmse/nsamps)
    
    bia_bypoint = bia_bypoint / float(nactdates)
    mae_bypoint = mae_bypoint / float(nactdates)
    rmse_bypoint = np.sqrt(rmse_bypoint/nactdates)
    
    ind = np.unravel_index(np.argmax(rmse_bypoint, axis=None), rmse_bypoint.shape)
    print ('   max rmse of ',rmse_bypoint[ind[0],ind[1]], ' at ',\
        ind[0], ind[1], lons[ind[0], ind[1]]-360., lats[ind[0], ind[1]])   
    #print ('nactdates = ', nactdates)
    
    # ---- calculate correlation coefficient and stddev for Taylor diagrams
    
    print ('   calculating corr coef and normalized stddev for Taylor diagrams. ')
    for jlat in range(nlats):
        for ilon in range(nlons):
            
            # ---- JFM 
            
            o = analysis_3d_prime[ibracket:91,jlat,ilon]
            f = forecast_3d_bcorr_prime[ibracket:91,jlat,ilon]
            r = np.corrcoef(o,f)
            corr_winter[jlat,ilon] = r[0,1]
            ostd = np.std(o)
            fstd = np.std(f)
            std_normalized_winter[jlat,ilon] = fstd/ostd
            std_obs_winter[jlat,ilon] = ostd
            
            # ---- AMJ 
            
            o = analysis_3d_prime[91:182,jlat,ilon]
            f = forecast_3d_bcorr_prime[91:182,jlat,ilon]
            r = np.corrcoef(o,f)
            corr_spring[jlat,ilon] = r[0,1]
            ostd = np.std(o)
            fstd = np.std(f)
            std_normalized_spring[jlat,ilon] = fstd/ostd
            std_obs_spring[jlat,ilon] = ostd
            
            # ---- JAS
            
            o = analysis_3d_prime[182:273,jlat,ilon]
            f = forecast_3d_bcorr_prime[182:273,jlat,ilon]
            r = np.corrcoef(o,f)
            corr_summer[jlat,ilon] = r[0,1]
            ostd = np.std(o)
            fstd = np.std(f)
            std_normalized_summer[jlat,ilon] = fstd/ostd
            std_obs_summer[jlat,ilon] = ostd
            
            # ---- OND 
            
            o = analysis_3d_prime[273:,jlat,ilon]
            f = forecast_3d_bcorr_prime[273:,jlat,ilon]
            r = np.corrcoef(o,f)
            corr_fall[jlat,ilon] = r[0,1]
            ostd = np.std(o)
            fstd = np.std(f)
            std_normalized_fall[jlat,ilon] = fstd/ostd
            std_obs_fall[jlat,ilon] = ostd
    
            #if std_normalized_summer[jlat,ilon] < 0.42 and corr_summer[jlat,ilon] < 0.1:
            #    n = len(o)
            #    for i in range(n):
            #        print ('i,o,f = ',i,o[i],f[i])
            #    print ('std corr = ', std_normalized_summer[jlat,ilon], corr_summer[jlat,ilon])
            #    print ('jlat,ilon = ', jlat, ilon, lons[jlat, ilon]-360., lats[jlat, ilon])
            #    sys.exit()
            
    # --- write to text file, cPickle file, and file for Taylor stuff
    
    ouf = open(statsfile, 'w')
    print (bia, mae, rmse, file=ouf)
    ouf.close()
    
    statsfile2 = statsfile[0:-4]+'.cPick'
    print ('   writing to ', statsfile2)
    ouf = open(statsfile2, 'wb')
    cPickle.dump(rmse_byday, ouf)
    cPickle.dump(bia_byday, ouf)
    cPickle.dump(mae_byday, ouf)
    cPickle.dump(rmse_bypoint, ouf)
    cPickle.dump(bia_bypoint, ouf)
    cPickle.dump(mae_bypoint, ouf)
    ouf.close()
    
    statsfile3 = statsfile[0:-4]+'_taylor.cPick'
    print ('   writing to ', statsfile3)
    ouf = open(statsfile3, 'wb')
    cPickle.dump(corr_winter, ouf)
    cPickle.dump(corr_spring, ouf)
    cPickle.dump(corr_summer, ouf)
    cPickle.dump(corr_fall, ouf)
    cPickle.dump(std_normalized_winter, ouf)
    cPickle.dump(std_normalized_spring, ouf)
    cPickle.dump(std_normalized_summer, ouf)
    cPickle.dump(std_normalized_fall, ouf)
    cPickle.dump(std_obs_winter,ouf)
    cPickle.dump(std_obs_spring,ouf)
    cPickle.dump(std_obs_summer,ouf)
    cPickle.dump(std_obs_fall,ouf)
    
    
    return rmse, bia, mae