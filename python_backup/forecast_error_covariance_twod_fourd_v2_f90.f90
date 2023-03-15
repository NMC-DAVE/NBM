SUBROUTINE forecast_error_covariance_twod_fourd_v2_f90( difference_threed, &
	efold, npts, ndates, nlats, nlons, Bx_localized, Bx_localized_fourd)
	
! sudo f2py --opt='-O4' -c -m forecast_error_covariance_twod_fourd_v2_f90 forecast_error_covariance_twod_fourd_v2_f90.f90 calculate_cov.f90
!
! forecast_error_covariance_twod_fourd_v2_f90.F90
! input is 3D array of bias-corrected forecasts (or bias corrections) across US. 
! Output is array of localized forecast-error covariances 
! (or bias-correction covariances), needed in Kalman filter.

INTEGER, INTENT(IN):: npts, ndates, nlats, nlons
REAL*8, INTENT(IN), DIMENSION(ndates, nlats, nlons) :: difference_threed
REAL*8, DIMENSION(npts, npts), INTENT(OUT) :: Bx_localized
REAL*8, DIMENSION(nlats, nlons, nlats, nlons), INTENT(OUT) :: Bx_localized_fourd

REAL*8 x1(ndates), x2(ndates)
REAL*8 cov
REAL*8 rmean
REAL*8 localizn_factor
REAL*8 hdist

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
		!difference_threed_prime(:,j1,i1) = difference_threed(:,j1,i1)
	END DO
END DO 
   
Bx_localized(:,:) = 0.0
Bx_localized_fourd(:,:,:,:) = 0.0

ktr1 = 1
DO j1 = 1, nlats
	DO i1 = 1, nlons
        x1(:) = difference_threed_prime(:,j1,i1)
        ktr2 = 1
       	DO j2 = 1, nlats
			DO i2 = 1, nlons
				x2(:) = difference_threed_prime(:,j2,i2)
	            hdist = (111./2.) * SQRT( REAL(i1-i2)**2 + REAL(j1-j2)**2)
                localizn_factor = exp(-(hdist/efold)**2.0)
				CALL calculate_cov(x1,x2,ndates,cov)
				!IF (cov*localizn_factor > 1000.) THEN
				!	PRINT *,'i1, j1, i2, j2 = ', i1, j1, i2, j2
				!	PRINT *,'sum(x1), sum(x2) = ', sum(x1), sum(x2)
				!	PRINT *,'ndates, cov, localizn_factor = ', ndates, cov, localizn_factor
				!END IF
                Bx_localized(ktr1,ktr2) = cov*localizn_factor
                Bx_localized(ktr2,ktr1) = Bx_localized(ktr1,ktr2)
                Bx_localized_fourD(j1,i1,j2,i2) = Bx_localized(ktr1,ktr2)
                Bx_localized_fourD(j2,i2,j1,i1) = Bx_localized(ktr1,ktr2)
                ktr2 = ktr2 + 1
			END DO
		END DO
        ktr1 = ktr1 + 1
	END DO
END DO
PRINT *, '   max, min B_localized = ', maxval(Bx_localized), minval(Bx_localized)

RETURN 
END SUBROUTINE forecast_error_covariance_twod_fourd_v2_f90
