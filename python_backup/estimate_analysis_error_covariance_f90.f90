SUBROUTINE estimate_analysis_error_covariance_f90(npts, nlats, nlons, &
	anal_err_var, efold_analysis, R)
	
! sudo f2py --opt='-O4' -c -m estimate_analysis_error_covariance_f90 estimate_analysis_error_covariance_f90.f90    
!
! estimate_analysis_error_covariance_f90.f90
!

INTEGER, INTENT(IN):: npts, nlons, nlats
REAL*8, INTENT(IN) :: anal_err_var, efold_analysis
REAL*8, DIMENSION(npts, npts), INTENT(OUT) :: R

REAL*8 localizn_factor 
REAL*8 cov

! f2py intent(in) npts, nlats, nlons
! f2py intent(in) anal_err_var, efold_analysis
! f2py depend(npts,npts) R
! f2py intent(out) R

ktr1 = 1
DO j1 = 1, nlats
	DO i1 = 1, nlons
        ktr2 = 1
		DO j2 = 1, nlats
        	DO i2 = 1, nlons
				hdist = (111./2.) * SQRT( REAL(i1-i2)**2 + REAL(j1-j2)**2)
                localizn_factor = exp(-(hdist/efold_analysis)**2.0)
                R(ktr1,ktr2) = anal_err_var*localizn_factor
                R(ktr2,ktr1) = anal_err_var*localizn_factor
                ktr2 = ktr2 + 1
			END DO
		END DO
        ktr1 = ktr1 + 1
	END DO
END DO

RETURN 
END SUBROUTINE estimate_analysis_error_covariance_f90
	