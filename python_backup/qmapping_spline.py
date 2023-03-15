def qmapping_spline(gefsv12_quantiles_on_ndfd, \
	precip_gefsv12_on_ndfd, spline_info_inv, \
	fraction_zero_ndfd, usegamma, 
	ny_ndfd, nx_ndfd, qmapped_precip):
	
import scipy.stats.gamma as gamma



REAL*8 alpha, beta, scale, cum, qgamma
REAL knots(17), bspline_coef(17), qmp, hazard_fn(1)

LOGICAL tootrue, toofalse

tootrue = True
toofalse = True

qmapped_precip = np.zeros((ny_ndfd,nx_ndfd), dtype=np.float32)
for jy in range(ny_ndfd):
    for ix in range(nx_ndfd):

	    if precip_gefsv12_on_ndfd[jy,ix] == 0.0:

			# ---- arbitrarily assign the CDF to zero if precip is zero.
        
			qmapped_precip[jy,ix] = 0.0
			 
		elif gefsv12_quantiles_on_ndfd[jy,ix] < fraction_zero_ndfd[jy,ix]:

			qmapped_precip(jy,ix)  = 0.0
            
		else:
            
			qpositive = (gefsv12_quantiles_on_ndfd[jy,ix] - fraction_zero_ndfd[jy,ix])/  \
		            (1.0 - fraction_zero_ndfd[jy,ix])

		    if usegamma[jy,ix] == 1:

		    	# ---- this was flagged as a dry point that estimated CDF with a Gamma.

		        alpha = spline_info_inv[jy,ix,0,0]   # previously stored here
		        beta = spline_info_inv[jy,ix,1,0]    # previously stored here
		        scale = 1.0 / beta
				cum = gefsv12_quantiles_on_ndfd[jy,ix]            
                qmp = gamma.ppf(qmp, alpha, loc=0, scale=scale)
				
				# --- with lack of training data, let's constrain the amount of quantile
				#     mapping possible
				
				if  qmp / precip_gefsv12_on_ndfd[jy,ix] .gt. 3.0 &
					qmp = precip_gefsv12_on_ndfd(jy,ix)*3.0
			ELSE
                
				IF (spline_info_inv(jy,ix,1,1) .eq. -99.99) THEN ! flag for no + training data
					
					qmp = precip_gefsv12_on_ndfd(jy,ix)

				ELSE
					
					! ---- flagged as a wet-enough point to estimate the CDF with 
					!      the spline fit to a hazard function. 
            
					knots(:) = spline_info_inv(jy,ix,1,:)
					bspline_coef(:) = spline_info_inv(jy,ix,2,:)
					hazard_fn(1) = -log(1.0 - qpositive)    
					CALL splev(knots, 17, bspline_coef, 3, hazard_fn, qmp, 1,ier)
					
				ENDIF
			ENDIF
			qmapped_precip(jy,ix) = qmp
		END IF
				
	END DO ! ix
END DO ! jy

RETURN
END SUBROUTINE qmapping_spline
