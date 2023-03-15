SUBROUTINE forecast_error_covariance_twod_fourd_f90( difference_threed, &
	npts, ndates, nlats, nlons, Bx_localized, Bx_localized_fourd)
	
! sudo f2py --opt='-O4' -c -m forecast_error_covariance_twod_fourd_f90 forecast_error_covariance_twod_fourd_f90.f90 calculate_var.f90
!
! forecast_error_covariance_twod_fourd_f90.F90
! input is 3D array of bias-corrected forecasts (or bias corrections) across US. 
! Output is array of localized forecast-error covariances 
! (or bias-correction covariances), needed in Kalman filter.

INTEGER, INTENT(IN):: npts, ndates, nlats, nlons
REAL*8, INTENT(IN), DIMENSION(ndates, nlats, nlons) :: difference_threed
REAL*8, DIMENSION(npts, npts), INTENT(OUT) :: Bx_localized
REAL*8, DIMENSION(nlats, nlons, nlats, nlons), INTENT(OUT) :: Bx_localized_fourd

REAL*8 x1(ndates)
REAL*8 cov
REAL*8 rmean

REAL*8, DIMENSION(ndates, nlats, nlons) :: difference_threed_prime

! f2py intent(in) npts, dates, nlats, nlons  # effsmall, efold, exponenty
! f2py intent(in) difference_threed
! f2py depend(ndates, nlats, nlons) difference_threed
! f2py intent(out) Bx_localized, Bx_localized_fourd
! f2py depend(npts, npts) Bx_localized
! f2py depend(nlats,nlons,nlats,nlons) Bx_localized_fourd
 
DO j1 = 1, nlats
	DO i1 = 1, nlons
        rmean = SUM(difference_threed(:,j1,i1)) / ndates
        difference_threed_prime(:,j1,i1) = difference_threed(:,j1,i1) - rmean
	END DO
END DO 
   
Bx_localized(:,:) = 0.0
Bx_localized_fourd(:,:,:,:) = 0.0
ktr1 = 1
DO j1 = 1, nlats
	DO i1 = 1, nlons
        x1(:) = difference_threed_prime(:,j1,i1)
		CALL calculate_var(x1,ndates,cov) 
        Bx_localized(ktr1,ktr1) = cov
        Bx_localized_fourd(j1,i1,j1,i1) = cov
        ktr1 = ktr1 + 1
	END DO
END DO

PRINT *, 'max, min Bx_localized = ', maxval(Bx_localized), minval(Bx_localized)

RETURN 
END SUBROUTINE forecast_error_covariance_twod_fourd_f90
