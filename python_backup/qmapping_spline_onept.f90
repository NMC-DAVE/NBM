SUBROUTINE qmapping_spline_onept(gefsv12_quantiles_on_ndfd, &
	precip_gefsv12_on_ndfd, spline_info_inv, &
	fraction_zero_ndfd, usegamma, qmapped_precip)
	
! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m qmapping_spline_onept qmapping_spline_onept.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f
	
REAL*8, INTENT(IN) :: gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, &
	fraction_zero_ndfd
INTEGER, INTENT(IN) :: usegamma
REAL*8, INTENT(IN), DIMENSION(2,17) :: spline_info_inv
REAL*8, INTENT(OUT) ::qmapped_precip

! f2py intent(in) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py intent(in) fraction_zero_ndfd, usegamma
! f2py intent(out) qmapped_precip
! f2py depend(2,17) spline_info_inv

REAL*8 alpha, beta, scale, cum, qgamma
REAL knots(17), bspline_coef(17), qmp, hazard_fn(1), qpositive

LOGICAL tootrue, toofalse

tootrue = .true.
toofalse = .false.

PRINT *,' gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, fraction_zero_ndfd, usegamma = ',&
	gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, fraction_zero_ndfd, usegamma

qmp= 0.0
IF (precip_gefsv12_on_ndfd .eq. 0.0) THEN

	! ---- arbitrarily assign the CDF to zero if precip is zero.

	qmapped_precip = 0.0
	 
ELSE IF (gefsv12_quantiles_on_ndfd .lt. fraction_zero_ndfd) THEN

	qmapped_precip  = 0.0
ELSE
	qpositive = (gefsv12_quantiles_on_ndfd - fraction_zero_ndfd)/  &
            (1.0 - fraction_zero_ndfd)
    IF (usegamma .eq. 1) THEN 

    	! ---- this was flagged as a dry point that estimated CDF with a Gamma.

        alpha = spline_info_inv(1,1)   ! previously stored here
        beta = spline_info_inv(2,1)    ! previously stored here
        scale = 1.0 / beta
		!cum = gefsv12_quantiles_on_ndfd
		
		! should above be qpositive
		!qmp = qgamma (cum,alpha,beta,tootrue,toofalse)
		qmp = qgamma (qpositive,alpha,beta,tootrue,toofalse)
		
		! --- with lack of training data, let's constrain the amount of quantile
		!     mapping possible
		
		IF (qmp/precip_gefsv12_on_ndfd .gt. 2.0) qmp=precip_gefsv12_on_ndfd*2.0
			
	ELSE
        
		IF (spline_info_inv(1,1) .eq. -99.99) THEN ! flag for no + training data
			
			qmp = precip_gefsv12_on_ndfd

		ELSE
			
			! ---- flagged as a wet-enough point to estimate the CDF with 
			!      the spline fit to a hazard function. 
    
			knots(:) = spline_info_inv(1,:)
			!print *,'knots = ', knots
			bspline_coef(:) = spline_info_inv(2,:)
			!print *,'bspline_coef = ', bspline_coef
			hazard_fn(1) = -log(1.0 - qpositive) 
			CALL splev(knots, 17, bspline_coef, 3, hazard_fn, qmp, 1,ier)
			PRINT *, 'qpositive, hazard_fn, qmp = ', qpositive, hazard_fn, qmp
			
		ENDIF
	ENDIF
	
	qmapped_precip = qmp
	
END IF
				

RETURN
END SUBROUTINE qmapping_spline_onept
