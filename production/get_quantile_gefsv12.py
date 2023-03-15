def get_quantile_gefsv12(precip_amount, precip_gefsv12_on_ndfd, \
    spline_info_gefsv12, ny_gefsv12, nx_gefsv12, \
    fraction_zero_gefsv12, fraction_zero_gefsv12_on_ndfd, \
    usegamma_gefsv12, quantile_99_ndfd, use99, \
    number_knots_gefsv12, lons_1d_realtime, \
    lats_1d_realtime_flipud, lons_ndfd, lats_ndfd, jnd, ind):    

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period.  This is
    the quantile in the overall distribution, including zeros. """

    import numpy as np
    from scipy.interpolate import LSQUnivariateSpline, splrep, splev
    import scipy.stats as stats
    from mpl_toolkits.basemap import Basemap, interp

    
    # ---- first loop thru the GEFSv12 grid points and determine the cumulative
    #      percentile of the forecast in the distribution of positive amounts.
    
    #print ('   precip_gefsv12_on_ndfd[jnd,ind] = ', precip_gefsv12_on_ndfd[jnd,ind])
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float64)
    #print ('   finding quantiles in + distribution')
    for jy in range(ny_gefsv12):
        for ix in range(nx_gefsv12):
            offset_out = 0.0
            if precip_amount[jy,ix] == 0.0:
        
                # ---- arbitrarily assign the CDF to zero if precip is zero.
        
                quantile = 0.0
                qpositive = 0.0
            else:	
        
                if usegamma_gefsv12[jy,ix] == 0: 
            
        	        # ---- flagged as a wet-enough point to estimate the CDF with 
                    #      the spline fit to a hazard function. 
                    
                    # Lesley, changed here.
            
                    nk = number_knots_gefsv12[jy,ix]
                    splines_tuple = (spline_info_gefsv12[jy,ix,0,0:nk], \
                        spline_info_gefsv12[jy,ix,1,0:nk], 3)
                    pmaxtrain = spline_info_gefsv12[jy,ix,0,nk-1] # max precip in the training data set
                    if precip_amount[jy,ix] < pmaxtrain:
                        spline_hazard = splev(precip_amount[jy,ix], splines_tuple)
                        qpositive = 1.0 - np.exp(-spline_hazard)
                    else: # beyond training data; 
                        qpositive = 0.9999
                    
                else:
            
                    if usegamma_gefsv12[jy,ix] == -1:
                        # --- flagged as basically no training data.
                        qpositive = 0.0
                    else:  # --- flagged as minimal training data 
                        #        so use Gamma distribution fit
                        alpha_hat = spline_info_gefsv12[jy,ix,0,0] 
                        beta_hat = spline_info_gefsv12[jy,ix,1,0] 
                        y0 = precip_amount[jy,ix] / beta_hat
                        qpositive = stats.gamma.cdf(y0, alpha_hat)
                
            gefsv12_quantiles[jy,ix] = qpositive
            
    # --- interpolate to the NDFD grid
    
    #print ('   interpolating quantiles to NDFD grid')
    gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles) 
    gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  
        
    # ---- if the quantile is above the 99th percentile and use99 == True, 
    #      then set the quantile to 0.99, and determine the offset
    
    #print ('   truncating + quantiles above 0.99 to 0.99')
    ninetynines = 0.99*np.ones((ny_ndfd, nx_ndfd), dtype=np.float64)
    gefsv12_quantiles_on_ndfd = np.where( gefsv12_quantiles_on_ndfd < 0.99, 
        gefsv12_quantiles_on_ndfd, ninetynines)

    # ---- for points with quantiles >= 0.99, set the offset

    #print ('   determining the offset for points with + quantiles >= 0.99')
    zeros = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    offset_on_ndfd = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    offset_on_ndfd = np.where( gefsv12_quantiles_on_ndfd >= 0.99, \
        precip_gefsv12_on_ndfd - quantile_99_ndfd, zeros)

    # ---- change the output quantile from the percentile of positive 
    #      values to the percentile in the distribution including zeros
    
    #print ('   changing output quantiles to include zeros.')
    gefsv12_quantiles_on_ndfd = fraction_zero_gefsv12_on_ndfd + \
        (1.0 - fraction_zero_gefsv12_on_ndfd)*gefsv12_quantiles_on_ndfd

    return gefsv12_quantiles_on_ndfd, offset_on_ndfd    
