SUBROUTINE upscale_precip_f90(apcp_anal,jboxmin, jboxmax, iboxmin, iboxmax, & 
	ny_half, nx_half, ny_ccpa, nx_ccpa, apcp_anal_upscaled )
	
	! compile for python call using
	! f2py --opt='-O4' -c -m upscale_precip_f90 upscale_precip_f90.f90 
	
	INTEGER, INTENT(IN) :: ny_half, nx_half, ny_ccpa, nx_ccpa
	INTEGER, INTENT(IN), DIMENSION(ny_half, nx_half) :: &
		jboxmin, jboxmax, iboxmin, iboxmax
	REAL, INTENT(IN), DIMENSION(ny_ccpa, nx_ccpa) :: apcp_anal
	REAL, INTENT(OUT), DIMENSION(ny_half, nx_half) :: apcp_anal_upscaled
	
	
	INTEGER npts
	!f2py intent(in) ny_half, nx_half, ny_ccpa, nx_ccpa
	!f2py intent(in) jboxmin, jboxmax, iboxmin, iboxmax, apcp_anal
	!f2py intent(out) apcp_anal_upscaled
	!f2py depend(ny_ccpa, nx_ccpa) apcp_anal
	!f2py depend(ny_half, nx_half) jboxmin, jboxmax, iboxmin, iboxmax
	!f2py depend(ny_half, nx_half) apcp_anal_upscaled
	
!	DO jy = 1, ny_half
!	DO jy = ny_half/2, ny_half/2
!		DO ix = 1, nx_half
			
	DO jy = 40,40
!	DO jy = ny_half/2, ny_half/2
		DO ix = 138,138
			
			jmin = jboxmin(jy,ix)+1  ! convert from python indices
			jmax = jboxmax(jy,ix)+1
			imin = iboxmin(jy,ix)+1
			imax = iboxmax(jy,ix)+1

			!print *,'jy, ix, jmin, jmax, imin, imax = ', jy, ix, jmin, jmax, imin, imax 
			!IF (jmin .lt. 1 .or. jmax .gt. ny_ccpa .or. imin .lt. 1 .or. imax .gt. nx_ccpa ) THEN
			!	print *,'box outside valid range: jy, ix, jmin, jmax, ny_ccpa, imin, imax, nx_ccpa = ',&
			!		jy, ix, jmin, jmax, ny_ccpa, imin, imax, nx_ccpa
			!ENDIF
			IF (jmin .ge. 1 .and. imin .ge. 1 .and. jmax .ge. 1 .and. imax .ge. 1) THEN
				npts = (jmax-jmin+1)*(imax-imin+1)
				!print *,'npts. (jmax-jmin+1), (imax-imin+1) = ', npts, (jmax-jmin+1), (imax-imin+1)
				apcp_anal_upscaled(jy,ix) = &
					SUM(apcp_anal(jmin:jmax,imin:imax)) / npts
				IF (apcp_anal_upscaled(jy,ix) .lt. -0.1) THEN
					apcp_anal_upscaled(jy,ix) = -99.99
				ENDIF
				IF (apcp_anal_upscaled(jy,ix) .ge. -0.1 .and. apcp_anal_upscaled(jy,ix) .le. 0.0) THEN
					apcp_anal_upscaled(jy,ix) = 0.0
				ENDIF
				!print *,'apcp_anal_upscaled(jy,ix) = ', apcp_anal_upscaled(jy,ix), maxval(apcp_anal(jmin:jmax,imin:imax))
			ELSE
				apcp_anal_upscaled(jy,ix) = -99.99
			END IF
		END DO
	END DO
	RETURN
	END SUBROUTINE upscale_precip_f90