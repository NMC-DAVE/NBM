SUBROUTINE qmapping_spline_fivexfive_flexiknot_test(gefsv12_quantiles_on_ndfd, &
	precip_gefsv12_on_ndfd, spline_info_inv, fraction_zero_gefsv12_ndfd, &
    fraction_zero_ndfd, usegamma, use99, offset, number_knots, nstride, &
    ny_ndfd, nx_ndfd, qmapped_precip)

! compile with the following python command
    
! f2py --opt='-O4' --opt='-Wno-tabs' -c -m qmapping_spline_fivexfive_flexiknot_test qmapping_spline_fivexfive_flexiknot_test.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f
	
! Quantile mapping using splines that predict the analyzed precipitation 
! given the cumulative hazard function.  Exceptions in the case of
! climatologically dry points (use gamma distributions for quantile mapping)
! or no ensemble mean precip (report back zero as quantile-mapped value).  

INTEGER, INTENT(IN) :: ny_ndfd, nx_ndfd, nstride
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: &
	gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, &
	fraction_zero_gefsv12_ndfd, fraction_zero_ndfd, offset
INTEGER, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: usegamma, number_knots
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd,2,17) :: spline_info_inv
LOGICAL, INTENT(IN) :: use99
REAL*8, INTENT(OUT), DIMENSION(25, ny_ndfd, nx_ndfd) ::	qmapped_precip

! f2py intent(in) ny_ndfd, nx_ndfd, nstride
! f2py intent(in) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py intent(in) fraction_zero_gefsv12_ndfd, usegamma, offset, use99
! f2py intent(in) number_knots, fraction_zero_ndfd
! f2py intent(out) qmapped_precip
! f2py depend(ny_ndfd, nx_ndfd) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py depend(ny_ndfd, nx_ndfd) fraction_zero_gefsv12_ndfd, usegamma, offset
! f2py depend(ny_ndfd, nx_ndfd) number_knots, fraction_zero_ndfd
! f2py depend(25,ny_ndfd, nx_ndfd) qmapped_precip
! f2py depend(ny_ndfd, nx_ndfd,2,17) spline_info_inv

INTEGER jyuse, ixuse
REAL*8 alpha, beta, scale, cum, qgamma 
REAL knots(17), bspline_coef(17)  ! the maximum value of knots, including 6 bounding.
   ! (3 lower, 3 upper).  Actual number in number_knots
REAL qmp, hazard_fn(1), qpositive

LOGICAL tootrue, toofalse

tootrue = .true.
toofalse = .false.

! ---- process every grid point in the NDFD domain.

DO jy = 1, ny_ndfd

	! ---- set the min and max j location for stencil loop.
	
	jymin = jy - 2*nstride
	jymax = jy + 2*nstride
	
	DO ix = 1, nx_ndfd
			
		! ---- set the min and max i location for stencil loop.
		ixmin = ix - 2*nstride
		ixmax = ix + 2*nstride
		istencil = 1

        ! ---- loop over a set of stencil points in the 
		!      surrounding the grid point of interest. Handle the 
        !      case where we are the edge of the domain too.
        
		DO jyf = jymin, jymax, nstride
			jyuse = jyf
			IF (jyf .lt. 1) jyuse = 1
			IF (jyf .gt. ny_ndfd) jyuse = ny_ndfd
			DO ixf = ixmin, ixmax, nstride
				ixuse = ixf
				IF (ixf .lt. 1) ixuse = 1
				IF (ixf .gt. nx_ndfd) ixuse = nx_ndfd

				IF (precip_gefsv12_on_ndfd(jyuse,ixuse) .eq. 0.0) THEN
                    ! ---- arbitrarily assign the CDF to zero if precip is zero.
					qmp = 0.0
				ELSE IF (gefsv12_quantiles_on_ndfd(jyuse,ixuse) .lt. &
                fraction_zero_gefsv12_ndfd(jy,ix)) THEN
					qmp  = 0.0
				ELSE IF (fraction_zero_gefsv12_ndfd(jy,ix) .gt. 0.99) THEN
					
					! handle the case where it's so dry that we are best off 
                    ! "quantile-mapped" value simply the raw value.
					qmp = precip_gefsv12_on_ndfd(jy,ix)

				ELSE
                    ! ---- offset being greater than zero is an indication
                    !      to assume the forecast quantile at the stencil point
                    !      is beyond the 99th percentile of wet values.  In this  
                    !      case, make sure the forecast percentile is properly set 

                    IF (offset(jyuse,ixuse) .gt. 0.0) THEN

						! ---- set the forecast quantile including zero values and 
						!      bound wet quantile to 0.99 
						
                    	qfcst_includezeros = fraction_zero_gefsv12_ndfd(jyuse,ixuse) + &
                    		(1.0 - fraction_zero_gefsv12_ndfd(jyuse,ixuse))*0.99
                    ELSE
                        !   ---- forecast is not exceptionally large.  The value 
						!        gefsv12_quantiles_on_ndfd passed in here will 
						!        include zero values
                        !        
                    	qfcst_includezeros = gefsv12_quantiles_on_ndfd(jyuse,ixuse)
                    END IF
                    
					! ---- for quantile mapping, the analyzed quantile is the same
					!      as the forecast quantile.
					
                    qanal_includezeros = qfcst_includezeros
					
					! ---- .... but the analyzed quantile of wet points (qpositive) 
					!      needs to back out the fraction zero of analyzed values.
					
                    qpositive = (qanal_includezeros - fraction_zero_ndfd(jy,ix)) /  &
                        (1.0 - fraction_zero_ndfd(jy,ix))

					IF (qpositive .lt. 0.0) qpositive = 0.0  ! just in case
							
	    			IF (usegamma(jy,ix) .eq. 1) THEN 

	    				! ---- this was flagged as a dry point that estimated 
						!      CDF with a Gamma distribution fit instead.
						
						alpha = spline_info_inv(jy,ix,1,1)   ! previously stored here
	        			beta = spline_info_inv(jy,ix,2,1)    ! previously stored here
						cum = qpositive ! gefsv12_quantiles_on_ndfd(jyuse,ixuse)
						qmp = qgamma(cum,alpha,beta,tootrue,toofalse)
			
						! --- with lack of training data, let's constrain the amount of quantile
						!     mapping possible to double that of the raw ensemble.
						
						IF (qmp / precip_gefsv12_on_ndfd(jyuse,ixuse) .gt. 2.0) &
							qmp = precip_gefsv12_on_ndfd(jyuse,ixuse)*2.0
					
					ELSE ! --- use splines
            
						IF (spline_info_inv(jy,ix,1,1) .eq. -99.99) THEN ! flag for problem
							qmp = precip_gefsv12_on_ndfd(jyuse,ixuse)
						ELSE
							! ---- flagged as a wet-enough point to estimate the 
							!      analyzed precipitation with the spline fit 
							!      predicting the precipitation amount given the cumulative
							!      Hazard function value.
                            
							knots(:) = spline_info_inv(jy,ix,1,:)
							bspline_coef(:) = spline_info_inv(jy,ix,2,:)
							hazard_fn(1) = -log(1.0 - qpositive) 
                            nk = number_knots(jy,ix)
							CALL splev(knots, nk, bspline_coef, 3, hazard_fn, qmp, 1,ier)
							
						ENDIF

					ENDIF
                    
					! ---- and finally, if today's precipitation amount was beyond the 
					!      99th percentile, we will add the previously determined offset
					!      which indicates how much greater the raw forecast was than the
					!      99th percentile.  This limits extreme quantile mappings in the
					!      tail of the distribution.
					
					IF (use99 .eqv. .TRUE. .and. offset(jy,ix) .gt. 0.0) THEN
						qmp = qmp + offset(jy,ix)
					ENDIF
                    
				ENDIF
				
				! ---- check for NaN's or below zero and fix.
				
				IF (isnan(qmp))	qmp = precip_gefsv12_on_ndfd(jyuse,ixuse)
				qmp=MAX(0.0, qmp)
				
				! ---- copy to output array.  With a 5x5 stencil of points, we want
				!      to keep track of each of these 25-fold greater ensembles of
				!      quantile mapped values. 
				
				qmapped_precip(istencil,jy,ix) = qmp
				istencil = istencil+1
				
			END DO ! ixf stencil point
		END DO ! jyf stencil point
				
	END DO ! ix grid point
END DO ! jy grid point


RETURN
END SUBROUTINE qmapping_spline_fivexfive_flexiknot_test
