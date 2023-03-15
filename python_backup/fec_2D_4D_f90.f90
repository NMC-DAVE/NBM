SUBROUTINE fec_2D_4D_f90a( difference_3d, efold, &
	exponenty, npts, ndates, nlats, nlons, Bx_localized, Bx_localized_4D)
	
! sudo f2py --opt='-O4' -c -m fec_2D_4D_f90 fec_2D_4D_f90.f90 
! sudo f2py -c -m fec_2D_4D_f90 fec_2D_4D_f90.f90       
!
! fec_2D_4D_f90.F90
! input is 3D array of bias-corrected forecasts (or bias corrections) across US. 
! Output is array of localized forecast-error covariances 
! (or bias-correction covariances), needed in Kalman filter.

INTEGER, INTENT(IN):: ndates, nlats, nlons, npts
REAL*8, INTENT(IN), DIMENSION(ndates, nlats, nlons) :: difference_3d
REAL*8, INTENT(IN) :: efold, exponenty
REAL*8, DIMENSION(npts, npts), INTENT(OUT) :: Bx_localized
REAL*8, DIMENSION(nlats, nlons, nlats, nlons), INTENT(OUT) :: Bx_localized_4D

REAL*8 x1(ndates)
REAL*8 x2(ndates)
REAL*8 localizn_factor 
REAL*8 cov

! f2py intent(in) ndates, nlats, nlons, npts
! f2py intent(in) difference_3d
! f2py intent(in) efold, exponenty
! f2py depend(ndates, nlats, nlons) difference_3d
! f2py intent(out) Bx_localized, Bx_localized_4D
! f2py depend(npts, npts) Bx_localized
! f2py depend(nlats,nlons,nlats,nlons) Bx_localized_4D
    
ktr1 = 1
DO i1 = 1, nlons
	DO j1 = 1, nlats
        x1(:) = difference_3d(:,j1,i1)
        ktr2 = 1
        DO i2 = 1, nlons
            DO j2 = 1, nlats
				x2(:) = difference_3d(:,j2,i2)
                !hdist = (111./2.) * sqrt( ( (i1-i2) * (coslat(j1)+coslat(j2)) / 2. )**2 + &
                !    (j1-j2)**2)
	            hdist = (111./2.) * SQRT( REAL(i1-i2)**2 + REAL(j1-j2)**2)
                localizn_factor = exp(-(hdist/efold)**exponenty)
				CALL calculate_cov(x1,x2,ndates,cov)
                Bx_localized(ktr1,ktr2) = cov*localizn_factor
                Bx_localized(ktr2,ktr1) = Bx_localized(ktr1,ktr2)
                Bx_localized_4D(j1,i1,j2,i2) = Bx_localized(ktr1,ktr2)
                Bx_localized_4D(j2,i2,j1,i1) = Bx_localized(ktr1,ktr2)
                ktr2 = ktr2 + 1
			END DO
		END DO
        ktr1 = ktr1 + 1
	END DO
END DO


RETURN 
END SUBROUTINE fec_2D_4D_f90a
	
! ======================================================================

SUBROUTINE calculate_cov(x1, x2, n, cov)
INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN), DIMENSION(n) :: x1, x2
REAL*8, INTENT(OUT) :: cov
REAL*8 x1mean, x2mean

x1mean = SUM(x1) / n
x2mean = SUM(x2) / n
cov = SUM( (x1(:)-x1mean) * (x2(:)-x2mean) ) / (n-1)
RETURN 
END SUBROUTINE calculate_cov

! ======================================================================

