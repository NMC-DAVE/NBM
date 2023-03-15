SUBROUTINE qmapping_spline_flexiknot_test(gefsv12_quantiles_on_ndfd, &
	precip_gefsv12_on_ndfd, spline_info_inv, fraction_zero_gefsv12_ndfd, &
    fraction_zero_ndfd, usegamma, use99, offset, number_knots, ny_ndfd, nx_ndfd, &
    qmapped_precip)

! compile with the following command
    
! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m qmapping_spline_flexiknot_test qmapping_spline_flexiknot_test.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f
	
! this version of the quantile mapping does a few things different from previous version.
! It is, per older articles such as Scheuerer and Hamill 20125, is going to use a stencil of 
! the surrounding grid points and their quantiles as expanded input into the quantile mapping.  

INTEGER, INTENT(IN) :: ny_ndfd, nx_ndfd
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: &
	gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, &
	fraction_zero_gefsv12_ndfd, fraction_zero_ndfd, offset
INTEGER, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: usegamma, number_knots
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd,2,17) :: spline_info_inv
LOGICAL, INTENT(IN) :: use99
REAL*8, INTENT(OUT), DIMENSION(ny_ndfd, nx_ndfd) ::	qmapped_precip

! f2py intent(in) ny_ndfd, nx_ndfd, nstride
! f2py intent(in) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py intent(in) fraction_zero_gefsv12_ndfd, fraction_zero_ndfd, 
! f2py intent(in) usegamma, offset, use99, number_knots 
! f2py intent(out) qmapped_precip
! f2py depend(ny_ndfd, nx_ndfd) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py depend(ny_ndfd, nx_ndfd) fraction_zero_gefsv12_ndfd, usegamma, offset
! f2py depend(ny_ndfd, nx_ndfd) number_knots, fraction_zero_ndfd
! f2py depend(ny_ndfd, nx_ndfd) qmapped_precip
! f2py depend(ny_ndfd, nx_ndfd,2,17) spline_info_inv

REAL*8 alpha, beta, scale, cum, qgamma 
REAL knots(17), bspline_coef(17)  ! the maximum value of knots.  Actual number in number_knots
REAL qmp, hazard_fn(1), qpositive, qfcst_includezeros, qfcst_positive, qanal_includezeros
INTEGER nz
LOGICAL tootrue, toofalse

tootrue = .true.
toofalse = .false.

DO jy = 1, ny_ndfd
	DO ix = 1, nx_ndfd
		
!DO jy = 480,480
!	DO ix = 696,696

		!PRINT *,'jy, ix = ',jy,ix
		IF (precip_gefsv12_on_ndfd(jy,ix) .eq. 0.0) THEN
			! ---- arbitrarily assign the CDF to zero if precip is zero.
			qmp = 0.0
		ELSE IF (gefsv12_quantiles_on_ndfd(jy,ix) .lt. &
		fraction_zero_gefsv12_ndfd(jy,ix)) THEN
			qmp  = 0.0
		ELSE IF (fraction_zero_gefsv12_ndfd(jy,ix) .gt. 0.99) THEN
			! handle the case where it's so dry that we are best off 
			! "quantile-mapped" value simply the raw value.
			qmp = precip_gefsv12_on_ndfd(jy,ix)
		ELSE
			! ---- offset being greater than zero is an indication 
			!      to assume the forecast quantile is beyond the 99th percentile.
			
			IF (offset(jy,ix) .gt. 0.0) THEN
				qfcst_positive = 0.99 
				qfcst_includezeros = fraction_zero_gefsv12_ndfd(jy,ix) + &
					(1.0 - fraction_zero_gefsv12_ndfd(jy,ix))*0.99 
			ELSE	
				qfcst_positive = (gefsv12_quantiles_on_ndfd(jy,ix) - &
					fraction_zero_gefsv12_ndfd(jy,ix)) /  &
					(1.0 - fraction_zero_gefsv12_ndfd(jy,ix))
				qfcst_includezeros = gefsv12_quantiles_on_ndfd(jy,ix)	
				
			END IF
			!IF (jy .eq. 480 .and. ix .eq. 696) &
			!	PRINT *,'qfcst_positive, qfcst_includezeros = ', qfcst_positive, qfcst_includezeros
			qanal_includezeros = qfcst_includezeros
			qpositive = (qanal_includezeros - fraction_zero_ndfd(jy,ix)) /  &
	            (1.0 - fraction_zero_ndfd(jy,ix))
			IF (qpositive .lt. 0.0) qpositive = 0.0	
				
			!IF (qpositive .lt. 0.0) THEN
			!	print *,'jy, ix = ', jy,ix
			!	print *,'qfcst_positive = ', qfcst_positive
			!	print *,'qfcst_includezeros = ', qfcst_includezeros
			!	print *,'qanal_includezeros = ', qanal_includezeros
			!	print *,'fraction_zero_gefsv12_ndfd(jy,ix) = ', fraction_zero_gefsv12_ndfd(jy,ix)
			!	print *,'fraction_zero_ndfd(jy,ix) = ', fraction_zero_ndfd(jy,ix)
			!	print *,'qpositive = ', qpositive
			!	print *,'alpha, beta, scale = ', alpha, beta, scale
			!	stop
			!ENDIF
					
			!IF (jy .eq. 480 .and. ix .eq. 696) &
			!	PRINT *, 'qanal_includezeros, qpositive = ', qanal_includezeros, qpositive
							
	    	IF (usegamma(jy,ix) .eq. 1) THEN 

	    		! ---- this was flagged as a dry point that estimated 
				!      CDF with a Gamma.
						
				alpha = spline_info_inv(jy,ix,1,1)   ! previously stored here
	        	beta = spline_info_inv(jy,ix,2,1)    ! previously stored here
	        	scale = 1.0 / beta
				cum = qpositive 
				qmp = qgamma(cum,alpha,beta,tootrue,toofalse)
			
				! --- with lack of training data, let's constrain the amount of quantile
				!     mapping possible
						
				IF (qmp / precip_gefsv12_on_ndfd(jy,ix) .gt. 2.0) &
					qmp = precip_gefsv12_on_ndfd(jy,ix)*2.0
			ELSE
            
				IF (spline_info_inv(jy,ix,1,1) .eq. -99.99) THEN ! flag for no + training data
					qmp = precip_gefsv12_on_ndfd(jy,ix)
				ELSE
					! ---- flagged as a wet-enough point to estimate the CDF with 
					!      the spline fit to a hazard function. 
                            
					knots(:) = spline_info_inv(jy,ix,1,:)
					bspline_coef(:) = spline_info_inv(jy,ix,2,:)
					hazard_fn(1) = -log(1.0 - qpositive) 
					IF (jy .eq. 480 .and. ix .eq. 696) PRINT *, '   hazard_fn, qpositive = ', hazard_fn, qpositive
					nk = number_knots(jy,ix)
					nz = 1
					CALL splev(knots, nk, bspline_coef, 3, hazard_fn, qmp, nz,ier)
					
		            !IF (jy .eq. 480 .and. ix .eq. 696) THEN       
					!	PRINT *,'   qmp after splev = ', qmp
					!ENDIF
                    
				ENDIF
				
			ENDIF
                    
			IF (use99 .eqv. .TRUE. .and. offset(jy,ix) .gt. 0.0) THEN
				qmp = qmp + offset(jy,ix)
			ENDIF
                    
            !IF (jy .eq. 480 .and. ix .eq. 696) THEN       
			!	print *,'   knots = ', knots(1:nk)
			!	print *,'   bspline_coef = ', bspline_coef(1:nk)
			!	PRINT *,'   qpositive = ', qpositive
			!	PRINT *,'   usegamma(459,603) = ', usegamma(480, 696)
			!	PRINT *,'   qmp after adding offset = ', qmp
			!ENDIF
                    
					
			!IF (qmp < -3.0) THEN
			!	PRINT *,'------ unexpected negative value at jy,ix = ', &
			!		jy,ix
			!	PRINT *,'ier, nk = ',ier, nk
			!	PRINT *,'precip_gefsv12_on_ndfd(jy,ix) =',&
			!		precip_gefsv12_on_ndfd(jy,ix)
			!	PRINT *,'qmp, usegamma, offset(jy,ix) = ', &
			!		qmp, usegamma(jy,ix), offset(jy,ix)
			!	print *,'gefsv12_quantiles_on_ndfd(jy,ix) = ',&
			!		gefsv12_quantiles_on_ndfd(jy,ix)
			!	print *, 'fraction_zero_gefsv12_ndfd(jy,ix) = ', &
			!		fraction_zero_gefsv12_ndfd(jy,ix)
			!	print *,'knots = ', knots(1:nk)
			!	print *,'bspline_coef = ', bspline_coef(1:nk)
			!	print *,'hazard_fn = ', hazard_fn
			!	print *,'qpositive = ', qpositive
			!	stop
			!ENDIF
		ENDIF
		
		IF (qmp < 0.0) qmp = 0.0		
		qmapped_precip(jy,ix) = qmp
				
	END DO ! ix
END DO ! jy

RETURN
END SUBROUTINE qmapping_spline_flexiknot_test
