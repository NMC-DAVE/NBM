SUBROUTINE update_beta(Kalman_gain_beta_4d, obsinc_2d, &
	idate, ndates, nlats, nlons, beta_3d)
	
	! sudo f2py --opt='-O4' -c -m update_beta update_beta.f90 
	
INTEGER, INTENT(IN) :: idate, ndates, nlats, nlons
REAL*8, INTENT(IN), DIMENSION(nlats, nlons, nlats, nlons) :: Kalman_gain_beta_4d
REAL*8, INTENT(IN), DIMENSION(nlats,nlons) :: obsinc_2d
REAL*8, INTENT(OUT), DIMENSION(ndates, nlats, nlons) :: beta_3d
	
! f2py intent(in) idate, ndates, nlats, nlons, Kalman_gain_beta_4d, obsinc_2d
! f2py depend(nlats, nlons, nlats, nlons)  Kalman_gain_beta_4d
! f2py depend(nlats, nlons) obsinc_2d
! f2py intent(out) beta_3d
! f2py depend(ndates, nlats, nlons) beta_3d

DO i = 1,nlons
    DO j = 1, nlats
        ktr = nlats*(j-1) + i
            
        ! ---- update the bias correction estimate, eq. 37 in Dee. 
            
        IF (idate > 1) THEN
            beta_3d(idate,j,i) = beta_3d(idate-1,j,i) - &
                SUM( Kalman_gain_beta_4D(j,i,:,:)*obsinc_2d(:,:) ) 
        ELSE
            beta_3d(idate,j,i) = - SUM (Kalman_gain_beta_4d(j,i,:,:)*obsinc_2d(:,:))
		END IF
	END DO
END DO
RETURN
END SUBROUTINE update_beta