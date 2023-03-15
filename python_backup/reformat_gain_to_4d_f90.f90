SUBROUTINE reformat_gain_to_4d_f90(gain_2D, nlats, nlons, gain_4D)

! sudo f2py --opt='-O4' -c -m reformat_gain_to_4d_f90 reformat_gain_to_4d_f90.f90 
	
INTEGER, INTENT(IN) :: nlats, nlons
REAL*8, INTENT(IN), DIMENSION(nlats*nlons, nlats*nlons) :: gain_2D
REAL*8, INTENT(OUT), DIMENSION(nlats,nlons,nlats,nlons) :: gain_4D
	
!f2py intent(in) gain_2D, nlats, nlons
!f2py intent(out) gain_4D
!f2py depend(nlats*nlons, nlats*nlons) :: gain_2D
!f2py depend(nlats,nlons,nlats,nlons) :: gain_4D

! ---- reform 2D Kalman gain into 4D-array.

ktr2 = 1
DO j2 = 1, nlats
	DO i2 = 1, nlons
		ktr1 = 1
		DO j1 = 1, nlats
			DO i1 = 1, nlons  
				!gain_4D(j1,i1,j2,i2) = gain_2D(ktr1,ktr2)
				gain_4D(j2,i2,j1,i1) = gain_2D(ktr2,ktr1)
				ktr1 = ktr1 + 1
			END DO
		END DO
		ktr2 = ktr2 + 1
	END DO
END DO
RETURN
END SUBROUTINE reformat_gain_to_4D_f90