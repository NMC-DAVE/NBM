SUBROUTINE forecast_error_covariance_twod_fourd_f90( difference_3d, &
	ndates, nlats, nlons, npts, Bx_localized, Bx_localized_4D)
	
! sudo f2py --opt='-O4' -c -m forecast_error_covariance_twod_fourd_f90 forecast_error_covariance_twod_fourd_f90.f90    
!
! forecast_error_covariance_twod_fourd_f90.F90
! input is 3D array of bias-corrected forecasts (or bias corrections) across US. 
! Output is array of localized forecast-error covariances 
! (or bias-correction covariances), needed in Kalman filter.

!INTEGER, INTENT(IN):: effsmall, ndates, nlats, nlons, npts 
INTEGER, INTENT(IN):: ndates, nlats, nlons, npts 
REAL*8, INTENT(IN), DIMENSION(ndates, nlats, nlons) :: difference_3d
!REAL*8, INTENT(IN) :: efold, exponenty
REAL*8, DIMENSION(npts, npts), INTENT(OUT) :: Bx_localized
REAL*8, DIMENSION(nlats, nlons, nlats, nlons), INTENT(OUT) :: Bx_localized_4D

REAL*8 x1(ndates)
!REAL*8 x2(ndates)
!REAL*8 localizn_factor 
REAL*8 cov
REAL*8 rmean

REAL*8, DIMENSION(ndates, nlats, nlons) :: difference_3d_prime

! f2py intent(in) ndates, nlats, nlons, npts  # effsmall, efold, exponenty
! f2py intent(in) difference_3d
! f2py depend(ndates, nlats, nlons) difference_3d
! f2py intent(out) Bx_localized, Bx_localized_4D
! f2py depend(npts, npts) Bx_localized
! f2py depend(nlats,nlons,nlats,nlons) Bx_localized_4D
 
DO j1 = 1, nlats
	DO i1 = 1, nlons
        rmean = SUM(difference_3d(:,j1,i1)) / ndates
        difference_3d_prime(:,j1,i1) = difference_3d(:,j1,i1) - rmean
	END DO
END DO 
 
!print *,'effsmall = ',effsmall   
!ktr1 = 1
!DO j1 = 1, nlats
!	DO i1 = 1, nlons
!        x1(:) = difference_3d(:,j1,i1)
!        ktr2 = 1
!       	DO j2 = 1, nlats
!			DO i2 = 1, nlons
!				x2(:) = difference_3d_prime(:,j2,i2)
!                !hdist = (111./2.) * sqrt( ( (i1-i2) * (coslat(j1)+coslat(j2)) / 2. )**2 + &
!                !    (j1-j2)**2)
!	            hdist = (111./2.) * SQRT( REAL(i1-i2)**2 + REAL(j1-j2)**2)
!                localizn_factor = exp(-(hdist/efold)**exponenty)
!				CALL calculate_cov(x1,x2,ndates,cov,effsmall)
!                Bx_localized(ktr1,ktr2) = cov*localizn_factor
!				!Bx_localized(ktr1,ktr2) = cov
!                Bx_localized(ktr2,ktr1) = Bx_localized(ktr1,ktr2)
!                Bx_localized_4D(j1,i1,j2,i2) = Bx_localized(ktr1,ktr2)
!                Bx_localized_4D(j2,i2,j1,i1) = Bx_localized(ktr1,ktr2)
!                ktr2 = ktr2 + 1
!			END DO
!		END DO
!        ktr1 = ktr1 + 1
!	END DO
!END DO



!print *,'effsmall = ',effsmall   
Bx_localized(:,:) = 0.0
Bx_localized_4D(:,:,:,:) = 0.09
ktr1 = 1
DO j1 = 1, nlats
	DO i1 = 1, nlons
        x1(:) = difference_3d_prime(:,j1,i1)
		CALL calculate_var(x1,ndates,cov) !effsmall
        Bx_localized(ktr1,ktr1) = cov
        Bx_localized_4D(j1,i1,j1,i1) = cov
        ktr1 = ktr1 + 1
	END DO
END DO

PRINT *, 'max, min Bx_localized = ', maxval(Bx_localized), minval(Bx_localized)

RETURN 
END SUBROUTINE forecast_error_covariance_twod_fourd_f90
	
! ======================================================================

SUBROUTINE calculate_var(x1, n, cov) ! effsmall
INTEGER, INTENT(IN) :: n  !, effsmall
REAL*8, INTENT(IN), DIMENSION(n) :: x1
REAL*8, INTENT(OUT) :: cov
REAL*8 denom

! kludge here to reduce effective sample size 0.3 factor in denom

denom = n-1
!IF (effsmall .eq. 1) THEN
!	denom = 0.3*(n-1)
!ELSE
!	denom = n-1
!ENDIF
!print *,'denom = ', denom
cov = SUM(x1(:)**2) / denom
!print *,'cov = ', cov
RETURN 
END SUBROUTINE calculate_var

! ======================================================================

