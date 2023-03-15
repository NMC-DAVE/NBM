def qmapping_spline_fivexfive_flexiknot(gefsv12_quantiles_on_ndfd, \
	precip_gefsv12_on_ndfd, spline_info_inv, fraction_zero_gefsv12_ndfd, \
    usegamma, use99, offset, number_knots, nstride, ny_ndfd, nx_ndfd):

    import scipy.interpolate.splev as splev
    import scipy.stats.gamma as gamma

    tootrue =  True
    toofalse = False

    # ---- process every grid point in the NDFD domain.

    for jy in range(ny_ndfd):
	    jymin = jy - 2*nstride
	    jymax = jy + 2*nstride
    
	    for ix range(nx_ndfd):
		    ixmin = ix - 2*nstride
		    ixmax = ix + 2*nstride
		    istencil = 0

            # ---- loop over a set of stencil points in the 
            #      surrounding the grid point of interest.
            #      Handle the case where we are the edge of the domain too.
        
		    for jyf in range(jymin, jymax, nstride):
			    jyuse = jyf
			    if jyf < 0: jyuse = 0
                if jyf > ny_ndfd-1 : jyuse = ny_ndfd-1
			    for ixf in range(ixmin, ixmax, nstride):
			        ixuse = ixf
			        if ixf < 0: ixuse = 0
				    if ixf > nx_ndfd-1 : ixuse = nx_ndfd-1
                
                    if precip_gefsv12_on_ndfd[jyuse,ixuse] == 0.0:
                        # ---- arbitrarily assign the CDF to zero if precip is zero.
                        qmp = 0.0
                    elif gefsv12_quantiles_on_ndfd[jyuse,ixuse] < \
                    fraction_zero_gefsv12_ndfd[jy,ix]:
                        # ---- overall fcst quantile less than fraction zero, so 
                        #      set precipitation to zero.
					    qmp  = 0.0
                    elif fraction_zero_gefsv12_ndfd[jy,ix] > 0.99:
					    # ---- handle the case where it's so climatologicall dry  
                        #      that we are best off setting the "quantile-mapped" 
                        #      value to simply the raw value.  We don't trust any
                        #      quantile mappings.
					    qmp = precip_gefsv12_on_ndfd[jy,ix]
				    else:
					    # ---- offset being greater than zero is an indication 
					    #      to assume the forecast quantile is beyond the 
                        #      99th percentile.
					    if offset[jyuse,ixuse] > 0.0:
						    qpositive = 0.99 
					    else:	
						    qpositive = (gefsv12_quantiles_on_ndfd[jyuse,ixuse] - \
							    fraction_zero_gefsv12_ndfd[jy,ix]) /  \
	            			    (1.0 - fraction_zero_gefsv12_ndfd[jy,ix])
						if usegamma[jy,ix] == 1: 

	    				    # ---- this was flagged as a dry point that estimated 
						    #      CDF with a Gamma.
						
						    alpha = spline_info_inv[jy,ix,0,0]   # previously stored here
                            beta = spline_info_inv[jy,ix,1,0]    # previously stored here
                            scale = 1.0 / beta
						    cum = qpositive # gefsv12_quantiles_on_ndfd(jyuse,ixuse)
                            qmp = gamma.ppf(qpositive, alpha, loc=0, scale=scale)
			
						    # --- with lack of training data, let's constrain the amount of quantile
						    #     mapping possible
						
						    if qmp / precip_gefsv12_on_ndfd[jyuse,ixuse] > 2.0:
							    qmp = precip_gefsv12_on_ndfd[jyuse,ixuse]*2.0
				
                        else:
            
						    if spline_info_inv[jy,ix,0,0] == -99.99: # flag for no + training data
							    qmp = precip_gefsv12_on_ndfd[jyuse,ixuse]
						    else:
                            
							    # ---- flagged as a wet-enough point to estimate the CDF with 
							    #      the spline fit to a hazard function. 
                            
                                nk = number_knots[jy,ix]
							    knots = spline_info_inv[jy,ix,0,0:nk]
                                bspline_coef = spline_info_inv[jy,ix,2,0:nk]
                                spline_tuple = (knots, bspline_coef, 3)
							    hazard_fn = -np.log(1.0 - qpositive)
							    qmp = splev(hazard_fn, spline_tuple)
                                
                    if use99 == True and offset[jy,ix] > 0.0: 
						qmp = qmp + offset[jy,ix]
				
				    qmapped_precip[istencil,jy,ix] = qmp
				    istencil = istencil+1
				

    return qmapped_precip


