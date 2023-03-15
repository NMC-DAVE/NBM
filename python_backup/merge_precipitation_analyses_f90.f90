SUBROUTINE merge_precipitation_analyses_f90(precip_ccpa, &
    precip_mswep, landmask, nlats_ndfd, nlons_ndfd, precip_final)
	
! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m merge_precipitation_analyses_f90 merge_precipitation_analyses_f90.f90 
! sudo f2py -c -m merge_precipitation_analyses_f90 merge_precipitation_analyses_f90.f90 

	
INTEGER, INTENT(IN) :: nlats_ndfd, nlons_ndfd
REAL*8, INTENT(IN), DIMENSION(nlats_ndfd, nlons_ndfd) :: precip_ccpa
REAL*4, INTENT(IN), DIMENSION(nlats_ndfd, nlons_ndfd) :: precip_mswep
INTEGER, INTENT(IN), DIMENSION(nlats_ndfd, nlons_ndfd) :: landmask
REAL*8, INTENT(OUT), DIMENSION(nlats_ndfd, nlons_ndfd) :: precip_final

!f2py intent(in) nlats_ndfd, nlons_ndfd
!f2py intent(in) precip_ccpa, precip_mswep, landmask
!f2py intent(out) precip_final
!f2py depend(ny,nx) precip_ccpa, precip_mswep, landmask, precip_final

PRINT *, 'Before loop'
DO jy = 1, nlats_ndfd
	PRINT *,'jy = ', jy
	DO ix = 1, nlons_ndfd
		IF (precip_ccpa(jy,ix) .ge. 0.0 .and. landmask(jy,ix) .eq. 1) THEN
			precip_final(jy,ix) = precip_ccpa(jy,ix)
		ELSE IF (precip_ccpa(jy,ix) .lt. 0.0 .and. landmask(jy,ix) .eq. 1) THEN
			precip_final(jy,ix) = precip_mswep(jy,ix)
		ELSE IF (landmask(jy,ix) .eq. 0) THEN
			precip_final(jy,ix) = precip_mswep(jy,ix)
		ELSE
			print *,'jy,ix,ccpa,mswep,landmask = ',jy, ix,&
				precip_ccpa(jy,ix), precip_mswep(jy,ix), precip_final(jy,ix)
		END IF
	END DO
END DO

RETURN
END SUBROUTINE merge_precipitation_analyses_f90