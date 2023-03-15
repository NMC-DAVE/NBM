def det_params(ndates_valid, nlats, nlons, data_input):

    """ try various power transformations to determine which power transformed normal distribution
        best fits the empirical CDF 
    """

    import numpy as np
    import scipy
    import sys

    power = np.zeros((nlats,nlons), dtype=np.float64)
    mean = np.zeros((nlats,nlons), dtype=np.float64)
    stddev = np.zeros((nlats,nlons), dtype=np.float64)
    
    empirical_quantiles = np.zeros((19,nlats, nlons), dtype=np.float64)
    fitted_quantiles = np.zeros((19), dtype=np.float64)
    Dn_statistic = np.zeros((19,nlats, nlons), dtype=np.float64)
    rq = 0.05 + 0.05*np.arange(19)
    iquse = (rq*ndates_valid).astype(int)

    # ---- determine the empirical quantiles from the sample
    
    print ('determining empirical quantiles')
    for ilat in range(nlats):
        for ilon in range(nlons):
            sample = data_input[:,ilat,ilon]
            sample_sorted = np.sort(sample)
            for iq in range(19):
                empirical_quantiles[iq,ilat,ilon] = sample_sorted[iquse[iq]]
            if ilat == 0 and ilon == 0: 
                print ('empirical quantiles = ', empirical_quantiles[:,ilat,ilon] )
                print ('sample_sorted= ', sample_sorted[0:-1:10])
                
    # ---- determine which power transformed normal distribution provides the best fit to 
    #      the empirical data
    
    print ('testing powers')
    testpowers = [0.5,0.7,0.8,0.9,0.95,1.0,1.05,1.1,1.2,1.3,1.5]
    ntest = len(testpowers)
    Dnstat = np.zeros((ntest), dtype=np.float64)
    Dnsample = np.zeros((19), dtype=np.float64)
    CDF_fitted = np.zeros((19), dtype=np.float64)
    for ilat in range(nlats):
        print ('ilat = ', ilat)
        for ilon in range(nlons):
            for itest, testpower in enumerate(testpowers):
                sample = empirical_quantiles[:,ilat,ilon]
                #print ('sample = ', sample)
                sample_xform = np.where(sample >= 0.0,\
                    ( (sample+1.0)**testpower - 1.0) / testpower, \
                    - ( (-sample+1.0)**(2.0-testpower) - 1.0 ) / (2.0-testpower) )
                smean = np.mean(sample_xform)
                sstd = np.std(sample_xform)
                #if ilat == 0 and ilon == 0: 
                #    print ('ilat,ilon, testpower, smean, sstd = ', ilat,ilon, testpower, smean, sstd)
                sample_norm = (sample_xform - smean) / sstd
                #print ('sample_norm = ', sample_norm)
                for iq in range(19):
                    CDF_fitted[iq] = scipy.stats.norm.cdf(sample_norm[iq], loc=0., scale=1.)
                    Dnsample[iq] = np.abs(CDF_fitted[iq] - rq[iq])
                #print ('Dnsample = ', Dnsample)
                Dnstat[itest] = np.max(Dnsample)
                    
            #print ('Dnstat = ', Dnstat)
            idmin = np.argmin(Dnstat)
            #print ('Dn min, index = ', Dnstat[idmin], idmin)
            testpower = testpowers[idmin]
            sample = data_input[:,ilat,ilon]
            sample_xform = np.where(sample >= 0.0,\
                ( (sample+1.0)**testpower - 1.0) / testpower, \
                - ( (-sample+1.0)**(2.0-testpower) - 1.0 ) / (2.0-testpower) )
            smean = np.mean(sample_xform)
            sstd = np.std(sample_xform)
            if ilat == 0 and ilon == 0: 
                print ('testpower, smean, sstd = ', testpower, smean, sstd)
            power[ilat,ilon] = testpower
            mean[ilat,ilon] = smean
            stddev[ilat,ilon] = sstd
            #sys.exit()
    print ('max, min power = ', np.max(power), np.min(power))
    print ('max, min mean = ', np.max(mean), np.min(mean))
    print ('max, min stddev = ', np.max(stddev), np.min(stddev))
    

    return power, mean, stddev