SUBROUTINE qmapping_spline(gefsv12_quantiles_on_ndfd, &
	precip_gefsv12_on_ndfd, spline_info_inv, &
	fraction_zero_ndfd, usegamma, use99, offset, &
	ny_ndfd, nx_ndfd, threshholds, nthresholds, nstride, &
	qmapped_precip, probability_forecast)
	
! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m qmapping_spline qmapping_spline.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f
	
INTEGER, INTENT(IN) ::  ny_ndfd, nx_ndfd, nthresholds, nstride
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: &
	gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, &
	fraction_zero_ndfd, offset
REAL*8, INTENT(IN), DIMENSION(nthresh) :: threshholds
INTEGER, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: usegamma
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd,2,17) :: spline_info_inv
LOGICAL, INTENT(IN) :: use99
REAL*8, INTENT(OUT), DIMENSION(ny_ndfd, nx_ndfd) ::	qmapped_precip
REAL*8, INTENT(OUT), DIMENSION(nthreshes, ny_ndfd, nx_ndfd) :: \
	probability_forecast

! f2py intent(in) ny_ndfd, nx_ndfd, nthresholds, nstride
! f2py intent(in) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py intent(in) fraction_zero_ndfd, usegamma, offset, use99
! f2py intent(in) thresholds
! f2py intent(out) qmapped_precip, probability_forecast
! f2py depend(ny_ndfd, nx_ndfd) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py depend(ny_ndfd, nx_ndfd) fraction_zero_ndfd, qmapped_precip, usegamma, offset
! f2py depend(ny_ndfd, nx_ndfd,2,17) spline_info_inv
! f2py depend(nthresh, ny_ndfd, nx_ndfd) probability_forecast

REAL*8 alpha, beta, scale, cum, qgamma, pnumer(nthresh), pdenom(nthresh)
REAL knots(17), bspline_coef(17), qmp, hazard_fn(1), qpositive

LOGICAL tootrue, toofalse

tootrue = .true.
toofalse = .false.

DO jy = 1, ny_ndfd
	
	jymin = jy - 2*nstride
	jymax = jy + 2*nstride
	
	DO ix = 1, nx_ndfd
		
		iymin = ix - 2*nstride
		iymax = ix + 2*nstride

		pnumer = 0.0
		pdenom = 0.0
		DO iyf = iymin, iymax
			IF (iyf .lt. 1) jyuse = 1
			IF (iyf .gt. ny_ndfd) jyuse = ny_ndfd
			DO ixf = jxmin, jxmax
				IF (ixf .lt. 1) ixuse = 1
				IF (ixf .gt. nx_ndfd) jxuse = nx_ndfd

				IF (precip_gefsv12_on_ndfd(jyuse,ixuse) .eq. 0.0) THEN

					! ---- arbitrarily assign the CDF to zero if precip is zero.
        
					qmp = 0.0
			 
				ELSE IF (gefsv12_quantiles_on_ndfd(jyuse,ixuse) .lt. fraction_zero_ndfd(jy,ix)) THEN

					qmp  = 0.0
			
				ELSE
			
					! ---- offset being greater than zero is an indication 
					!      to assume the quantile is beyond the 99th percentile. 
			
					IF (offset(jyuse,ixuse) .gt. 0.0) THEN
						qpositive = (0.99 - fraction_zero_ndfd(jy,ix))/  &
	            			(1.0 - fraction_zero_ndfd(jy,ix))
					ELSE	
						qpositive = (gefsv12_quantiles_on_ndfd(jyuse,ixuse) - &
							fraction_zero_ndfd(jy,ix)) /  &
	            			(1.0 - fraction_zero_ndfd(jy,ix))
					END IF
		
	    			IF (usegamma(jy,ix) .eq. 1) THEN 

	    				! ---- this was flagged as a dry point that estimated 
						!      CDF with a Gamma.

	        			alpha = spline_info_inv(jy,ix,1,1)   ! previously stored here
	        			beta = spline_info_inv(jy,ix,2,1)    ! previously stored here
	        			scale = 1.0 / beta
						cum = gefsv12_quantiles_on_ndfd(jyuse,ixuse)
						qmp = qgamma(cum,alpha,beta,tootrue,toofalse)
			
						! --- with lack of training data, let's constrain the amount of quantile
						!     mapping possible
			
						IF (qmp / precip_gefsv12_on_ndfd(jyuse,ixuse) .gt. 2.0) &
							qmp = precip_gefsv12_on_ndfd(jyuse,ixuse)*2.0
				
						ELSE
            
							IF (spline_info_inv(jy,ix,1,1) .eq. -99.99) THEN ! flag for no + training data
				
								qmp = precip_gefsv12_on_ndfd(jyuse,ixuse)

							ELSE
				
								! ---- flagged as a wet-enough point to estimate the CDF with 
								!      the spline fit to a hazard function. 
        
								knots(:) = spline_info_inv(jy,ix,1,:)
								bspline_coef(:) = spline_info_inv(jy,ix,2,:)
								hazard_fn(1) = -log(1.0 - qpositive) 
								CALL splev(knots, 17, bspline_coef, 3, hazard_fn, qmp, 1,ier)
				
							ENDIF
						ENDIF
					
					END IF
	
	
					IF (use99 .eqv. .TRUE. .and. offset(jyuse,ixuse) .gt. 0.0) THEN
						qmp = qmp + offset(jyuse,ixuse)
					ENDIF
					
					IF (jyuse .eq. jy .and. ixuse .eq. ix) qmapped_precip(jy,ix) = qmp
					
				ENDIF
				
				! --- increment probability counters
				
				DO ithresh = 1, nthresholds
					pdenom(ithresh) = pdenom(ithresh)+1
					IF (qmp .gt. thresholds(ithresh)) pnumer(ithresh) = pnumer(ithresh) + 1.0
				END DO
		
			END DO ! ixf
		END DO ! jyf

		! ---- make final probability for this member.
		
		DO ithresh = 1, nthresholds
			probability_forecast(ithresh,jy,ix) = pnumer(ithresh) / pdenom(ithresh)
		END DO
				
	END DO ! ix
END DO ! jy

RETURN
END SUBROUTINE qmapping_spline
