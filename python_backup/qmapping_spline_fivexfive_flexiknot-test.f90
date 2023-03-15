SUBROUTINE qmapping_spline_fivexfive_flexiknot(gefsv12_quantiles_on_ndfd, &
	precip_gefsv12_on_ndfd, spline_info_inv, fraction_zero_gefsv12_ndfd, &
    usegamma, use99, offset, number_knots, nstride, ny_ndfd, nx_ndfd, &
    qmapped_precip)

! compile with the following command
    
! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m qmapping_spline_fivexfive_flexiknot qmapping_spline_fivexfive_flexiknot.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f
	
! this version of the quantile mapping does a few things different from previous version.
! It is, per older articles such as Scheuerer and Hamill 20125, is going to use a stencil of 
! the surrounding grid points and their quantiles as expanded input into the quantile mapping.  

INTEGER, INTENT(IN) :: ny_ndfd, nx_ndfd, nstride
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: &
	gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd, &
	fraction_zero_gefsv12_ndfd, offset
INTEGER, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd) :: usegamma, number_knots
REAL*8, INTENT(IN), DIMENSION(ny_ndfd, nx_ndfd,2,17) :: spline_info_inv
LOGICAL, INTENT(IN) :: use99
REAL*8, INTENT(OUT), DIMENSION(25, ny_ndfd, nx_ndfd) ::	qmapped_precip

! f2py intent(in) ny_ndfd, nx_ndfd, nstride
! f2py intent(in) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py intent(in) fraction_zero_gefsv12_ndfd, usegamma, offset, use99
! f2py intent(in) number_knots
! f2py intent(out) qmapped_precip
! f2py depend(ny_ndfd, nx_ndfd) gefsv12_quantiles_on_ndfd, precip_gefsv12_on_ndfd
! f2py depend(ny_ndfd, nx_ndfd) fraction_zero_gefsv12_ndfd, usegamma, offset
! f2py depend(ny_ndfd, nx_ndfd) number_knots
! f2py depend(25,ny_ndfd, nx_ndfd) qmapped_precip
! f2py depend(ny_ndfd, nx_ndfd,2,17) spline_info_inv

INTEGER jyuse, ixuse
REAL*8 alpha, beta, scale, cum, qgamma 
REAL knots(17), bspline_coef(17)  ! the maximum value of knots.  Actual number in number_knots
REAL qmp, hazard_fn(1), qpositive

LOGICAL tootrue, toofalse

tootrue = .true.
toofalse = .false.

! process every grid point in the NDFD domain.

DO jy = 1, ny_ndfd
	jymin = jy - 2*nstride
	jymax = jy + 2*nstride
	
	DO ix = 1, nx_ndfd
			
		ixmin = ix - 2*nstride
		ixmax = ix + 2*nstride
		istencil = 1

        ! ---- loop over a set of stencil points in the surrounding the grid point of interest.
        !      Handle the case where we are the edge of the domain too.
        
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
                    !qmp = precip_gefsv12_on_ndfd(jyuse,ixuse)
				ELSE
					! ---- offset being greater than zero is an indication 
					!      to assume the forecast quantile is beyond the 99th percentile.
					IF (offset(jyuse,ixuse) .gt. 0.0) THEN
						qpositive = 0.99 
					ELSE	
						qpositive = (gefsv12_quantiles_on_ndfd(jyuse,ixuse) - &
							fraction_zero_gefsv12_ndfd(jy,ix)) /  &
	            			(1.0 - fraction_zero_gefsv12_ndfd(jy,ix))
					END IF
							
	    			IF (usegamma(jy,ix) .eq. 1) THEN 

	    				! ---- this was flagged as a dry point that estimated 
						!      CDF with a Gamma.
						
						alpha = spline_info_inv(jy,ix,1,1)   ! previously stored here
	        			beta = spline_info_inv(jy,ix,2,1)    ! previously stored here
	        			scale = 1.0 / beta
						cum = qpositive ! gefsv12_quantiles_on_ndfd(jyuse,ixuse)
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
                            nk = number_knots(jy,ix)
							CALL splev(knots, nk, bspline_coef, 3, hazard_fn, qmp, 1,ier)
							qmp=MAX(0.0, qmp)
						ENDIF
					ENDIF
                    
					IF (use99 .eqv. .TRUE. .and. offset(jy,ix) .gt. 0.0) THEN
						qmp = qmp + offset(jy,ix)
					ENDIF
                    
                    ! this is a bunch of junk diagnostic code for debugging.
                    
                    IF (qmp < -3.0) THEN
                        PRINT *,'------ unexpected negative value at jy,ix, jyuse, ixuse = ', &
                            jy,ix, jyuse, ixuse
                        PRINT *,'ier, nk = ',ier, nk
						PRINT *,'precip_gefsv12_on_ndfd(jyuse,ixuse) =',&
							precip_gefsv12_on_ndfd(jyuse,ixuse)
                        PRINT *,'qmp, usegamma, offset(jyuse,ixuse) = ', &
                            qmp, usegamma(jy,ix), offset(jyuse,ixuse)
                        print *,'gefsv12_quantiles_on_ndfd(jy,ix), (jyuse,ixuse) = ',&
                            gefsv12_quantiles_on_ndfd(jy,ix), &
                            gefsv12_quantiles_on_ndfd(jyuse,ixuse)
                        print *, 'fraction_zero_gefsv12_ndfd(jy,ix) = ', &
                            fraction_zero_gefsv12_ndfd(jy,ix)
                        print *,'knots = ', knots(1:nk)
                        print *,'bspline_coef = ', bspline_coef(1:nk)
                        print *,'hazard_fn = ', hazard_fn
                        print *,'qpositive = ', qpositive		
                        
                        DO iq = 1,99
                            rq = REAL(iq) / 100.
                            hazard_fn(1) = -log(1.0 - rq) 
                            CALL splev(knots, nk, bspline_coef, 3, hazard_fn, qmp, 1,ier)
                            PRINT *, 'rq, qmp = ', rq, qmp
                        END DO                        		
                        stop
                    ENDIF 
                    
				ENDIF
				
				qmapped_precip(istencil,jy,ix) = qmp
				istencil = istencil+1
				
			END DO ! ixf
		END DO ! jyf
				
	END DO ! ix
END DO ! jy


RETURN
END SUBROUTINE qmapping_spline_fivexfive_flexiknot
