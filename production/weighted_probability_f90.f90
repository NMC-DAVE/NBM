SUBROUTINE weighted_probability_f90(precip_mean, precip_ens, thresholds, &
	thresholds_hist, closest_hist, ny, nx, nmembers, nthresholds, nchcats, prob)

	! f2py --opt='-O4' -c -m weighted_probability_f90 weighted_probability_f90.f90
	! then cp ./weighted_probability_f90.cpython-38-darwin.so ./weighted_probability_f90.so 
	
	INTEGER, INTENT(IN) :: ny, nx, nmembers, nthresholds, nchcats
	REAL, INTENT(IN), DIMENSION(ny, nx) :: precip_mean
	REAL, INTENT(IN), DIMENSION(nthresholds) :: thresholds
	REAL, INTENT(IN), DIMENSION(nmembers, ny, nx) :: precip_ens
	REAL, INTENT(IN), DIMENSION(nchcats-1) :: thresholds_hist
	REAL*8, INTENT(IN), DIMENSION(nmembers, nchcats) :: closest_hist
	REAL*8, INTENT(OUT), DIMENSION(nthresholds, ny, nx) :: prob

	!f2py intent(in) ny, nx, nmembers, nthresholds, nchcats
	!f2py intent(in) precip_ens
	!f2py intent(in) thresholds, thresholds_hist
	!f2py intent(in) closest_hist
	!f2py intent(out) prob
	!f2py depend(ny,nx) precip_mean
	!f2py depend(nthresholds, ny,nx) prob
	!f2py depend(nmembers, ny, nx) precip_ens
	!f2py depend(nchcats-1) thresholds_hist
	!f2py depend(nthresholds) thresholds
	!f2py depend(nmembers, nchcats) closest_hist

	! --- local variables
	
	REAL*8 ens1d(nmembers), ens1d_sorted(nmembers), p
	INTEGER bines(nmembers)
	REAL t(nchcats-1)
	INTEGER maxlocs(2)
	
	! ---- loop over grid points
	
	!PRINT *,'pthresh = ', pthresh
	!PRINT *,'closest_hist[1:nmembers:10,1] = ', closest_hist(1:nmembers:10,1) 
	!PRINT *,'closest_hist[1:nmembers:10,2] = ', closest_hist(1:nmembers:10,2) 
	!PRINT *,'closest_hist[1:nmembers:10,3] = ', closest_hist(1:nmembers:10,3)
	!PRINT *,'closest_hist[1:nmembers:10,4] = ', closest_hist(1:nmembers:10,4)
	!PRINT *,'max, min precip_ens = ', MAXVAL(precip_ens), MINVAL(precip_ens)
	
	!print *,'ny, nx, nmembers, nthresholds, nchcats = ', ny, nx, nmembers, nthresholds, nchcats
	
	prob(:,:,:) = 0.0
	t = thresholds_hist
	
	pmax = 0.0
	DO jy = 1, ny
		DO ix = 1, nx
			IF (precip_mean(jy,ix) .gt. pmax) THEN
				pmax = precip_mean(jy,ix)
				imax = ix
				jmax = jy
			END IF
		END DO
	END DO

	PRINT *,'  pmax at jy, ix = ',precip_mean(jmax,imax), jmax,imax
	
	DO jy = 1, ny
		!IF (mod(jy,100) .eq. 0) PRINT *,'      jy = ',jy,' of ',ny
		DO ix = 1, nx
	!DO jy = jmax, jmax
	!	DO ix = imax -200, imax+200, 20
			
			! ---- determine the category for the closest-member
			!      histogram CDF to apply for this point based
			!      on quantile-mapped ensemble-mean precipitation
			
			
			p = precip_mean(jy,ix)
			IF (p .eq. 0.0) THEN
				prob(:,jy,ix) = 0.0 ! don't bother sorting and weighting
			ELSE
				IF (p .le. t(1)) THEN
					icat = 1 ! near zero precip
				ELSE IF (p .gt. t(1) .and. p .lt. t(2)) THEN
					icat = 2 ! light precip
				ELSE IF (p .gt. t(2) .and. p .lt. t(3)) THEN
					icat = 3 ! light-moderate precip
				ELSE IF (p .gt. t(3) .and. p .lt. t(4)) THEN
					icat = 4 ! moderate precip
				ELSE IF (p .gt. t(4) .and. p .lt. t(5)) THEN
					icat = 5 ! moderate-heavy precip
				ELSE
					icat = 6 ! heavy precip
				ENDIF
			
				! ---- count up the number of sorted ensemble members greater than
				!      or = to this threshold.  Sum the closest-member histogram 
				!      associated with those sorted members.
			
				ens1d(:) = precip_ens(:,jy,ix)
				DO ithresh = 1, nthresholds
					numgt = 0
					DO imem = 1, nmembers
						IF (ens1d(imem) .ge. thresholds(ithresh)) numgt = numgt+1
					END DO
					prob(ithresh,jy,ix) = SUM(closest_hist(nmembers-numgt+1:nmembers, icat))
				END DO
				
				!IF (jy .ge. jmax-20 .and. jy .le. jmax+20 .and. &
				!ix .ge. imax-20 .and. ix .le. imax+20) THEN
				!	PRINT *,'====== jy, ix, precip_mean = ', jy, ix, precip_mean(jy,ix)
				!	!print *,'    ens1d = ', ens1d
				!	print *,'    pthresh, numgt, nmembers = ', pthresh, numgt, nmembers
				!	!print *,'    closest_hist(nmembers-imem+1:nmembers, icat)'
				!	print *,'    SUM(closest_hist(nmembers-numgt+1:nmembers, icat))',&
				!		SUM(closest_hist(nmembers-numgt+1:nmembers, icat))
				!	print *,'    prob(jy,ix) = ', prob(jy,ix)
				!END IF
				
			END IF
			
			!PRINT *,'------ jy, ix = ',jy, ix
			!PRINT *,'   pmean, icat, numgt = ', precip_mean(jy,ix), icat, numgt
			!PRINT *,'   prob(jy,ix) = ', prob(jy,ix)
			!PRINT *,'   closest_hist(nmembers-numgt+1:nmembers, icat) = ', &
			!	closest_hist(nmembers-numgt+1:nmembers,icat )

		END DO
	END DO
	PRINT *,'  max, min prob = ', maxval(prob), minval(prob)
	!PRINT *,'prob(jmax,imax) = ', prob(jmax,imax) 
	!print *,'prob(jmax-10:jmax+10,imax) = ',prob(jmax-10:jmax+10,imax)
	!print *,'prob(jmax,imax-10:imax+10) = ',prob(jmax,imax-10:imax+10)

	RETURN
	
	END SUBROUTINE weighted_probability_f90


