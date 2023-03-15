SUBROUTINE reformat_2d_to_4d_f90(arr_2D, nlats, nlons, arr_4D)

! sudo f2py --opt='-O4' -c -m reformat_4d_to_2d_f90 reformat_4d_to_2d_f90.f90 
	
INTEGER, INTENT(IN) :: nlats, nlons
REAL*8, INTENT(OUT), DIMENSION(nlats,nlons,nlats,nlons) :: arr_4D
REAL*8, INTENT(IN), DIMENSION(nlats*nlons, nlats*nlons)  :: arr_2D
	
!f2py intent(in) arr_2D, nlats, nlons
!f2py intent(out) arr_4D
!f2py depend(nlats*nlons, nlats*nlons) :: arr_2D
!f2py depend(nlats,nlons,nlats,nlons) :: arr_4D

! ---- reform 4D array into 2D-array.

ktr1 = 1
DO i1 = 1, nlons  
	DO j1 = 1, nlats
		ktr2 = 1
		DO i2 = 1, nlons
			DO j2 = 1, nlats
				arr_4D(j1,i1,j2,i2) = arr_2D(ktr1,ktr2)
				ktr2 = ktr2 + 1
			END DO
		END DO
		ktr1 = ktr1 + 1
	END DO
END DO

RETURN
END SUBROUTINE reformat_2d_to_4d_f90