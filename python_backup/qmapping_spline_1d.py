def qmapping_spline_1d(gefsv12_quantiles_on_ndfd, \
	precip_gefsv12_on_ndfd, spline_info_inv, \
	fraction_zero_ndfd, usegamma):
	
    import scipy.stats as stats
    from scipy.interpolate import LSQUnivariateSpline, splrep, splev
    import numpy as np

    qmapped_precip = 0.0

    if precip_gefsv12_on_ndfd == 0.0:

    	# ---- arbitrarily assign the CDF to zero if precip is zero.

    	qmapped_precip = 0.0
	 
    elif gefsv12_quantiles_on_ndfd < fraction_zero_ndfd:

    	qmapped_precip  = 0.0
    
    else:
    
        qpositive = (gefsv12_quantiles_on_ndfd - fraction_zero_ndfd) /  \
            (1.0 - fraction_zero_ndfd)
        if usegamma == 1:

            # ---- this was flagged as a dry point that estimated CDF 
            #      with a Gamma.

            alpha = spline_info_inv[0,0]   # previously stored here
            beta = spline_info_inv[1,0]    # previously stored here
            scale = 1.0 / beta
            cum = gefsv12_quantiles_on_ndfd            
            qmp = stats.gamma.ppf(qpositive, alpha, loc=0, scale=scale)
            
            # --- with lack of training data, let's constrain the amount 
            #     of quantile mapping possible
            
            if  qmp / precip_gefsv12_on_ndfd > 2.0:
                qmp = precip_gefsv12_on_ndfd*2.0
                
        else:
        
            if spline_info_inv[0,0] < -99.9: # flag for no + training data
                qmp = precip_gefsv12_on_ndfd
            else:
                # ---- flagged as a wet-enough point to estimate the CDF with 
                #      the spline fit to a hazard function. 
                
                splines_tuple = (spline_info_inv[0,:], \
                    spline_info_inv[1,:], 3)
                hazard_fn = -np.log(1.0 - qpositive) 
                qmapped_precip = splev(hazard_fn, splines_tuple, ext=0)
                #print ('   qpositive, hazard_fn, qmapped_precip = ', \
                #    qpositive, hazard_fn, qmapped_precip)

    return qmapped_precip

