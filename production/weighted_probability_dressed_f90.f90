SUBROUTINE weighted_probability_dressed_f90(precip_mean, precip_ens, &
    thresholds, thresholds_hist, closest_hist, dress_array, &
	ny, nx, nmembers, nthresholds, ndress, nchcats, prob)
	
	! f2py --opt='-O2' -c -m weighted_probability_dressed_f90 weighted_probability_dressed_f90.f90 cumnor.f90 quicksort.f90 r8_swap.f90
	
	! then $ cp weighted_probability_dressed_f90.cpython-38-darwin.so weighted_probability_dressed_f90.so

	!integer, parameter :: sp = selected_real_kind(6,37) ! single precision
	!integer, parameter :: dp = selected_real_kind(15,307) ! double precision
	
	INTEGER*4, INTENT(IN) :: ny, nx, nmembers, nthresholds, nchcats, ndress
	REAL, INTENT(IN), DIMENSION(ny, nx) :: precip_mean
	REAL, INTENT(IN), DIMENSION(nmembers, ny, nx) :: precip_ens
	REAL, INTENT(IN), DIMENSION(nthresholds) :: thresholds
	REAL, INTENT(IN), DIMENSION(nchcats-1) :: thresholds_hist
	REAL, INTENT(IN), DIMENSION(ndress) :: dress_array

	REAL*8, INTENT(IN), DIMENSION(nmembers, nchcats) :: closest_hist
	REAL*8, INTENT(OUT), DIMENSION(nthresholds, ny, nx) :: prob

	!f2py intent(in) ny, nx, nmembers, nthresholds, ndress, nchcats
	!f2py intent(in) precip_ens, precip_mean 
	!f2py intent(in) thresholds, thresholds_hist
	!f2py intent(in) closest_hist
	!f2py intent(in) dress_array
	!f2py depend(nmembers, ny, nx) precip_ens
	!f2py depend(ny,nx) precip_mean
	!f2py depend(nthresholds) thresholds
	!f2py depend(nchcats-1) thresholds_hist
	!f2py depend(nmembers, nchcats) closest_hist
	!f2py depend(ndress) dress_array
	!f2py depend(nthresholds, ny,nx) prob
	!f2py intent(out) prob

	! --- local variables
	
	REAL*8 ens1d(nmembers), p, zscore, cum, ccum, pmean, pstd
	REAL t(nchcats-1)
	REAL*8 cum_table(-5000:5000)
	REAL b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, b1_std_lowrank 
	REAL b0_mean_midrank, b1_mean_midrank, b0_std_midrank, b1_std_midrank 
	REAL b0_mean_highrank, b1_mean_highrank, b0_std_highrank, b1_std_highrank 
	
	
	b0_mean_lowrank = dress_array(1)
	b1_mean_lowrank = dress_array(2)
	b0_std_lowrank = dress_array(3)
	b1_std_lowrank = dress_array(4)
	b0_mean_midrank = dress_array(5) 
	b1_mean_midrank = dress_array(6)
	b0_std_midrank = dress_array(7) 
	b1_std_midrank = dress_array(8)
	b0_mean_highrank = dress_array(9)
	b1_mean_highrank = dress_array(10)
	b0_std_highrank = dress_array(11)
	b1_std_highrank = dress_array(12)
	
	
	! ---- make a lookup table of cumulative exceedance probabilities
	
	DO i = -5000, 5000
		zscore = DBLE(i)/1000.
		CALL cumnor(zscore, cum, ccum) 
		cum_table(i) = cum
	END DO
			
	! ---- determine point with max ensemble mean precip	
	
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
	PRINT *,'  precip mean max at jy, ix = ',precip_mean(jmax,imax), jmax,imax
	
	! ---- process each grid point
	
	prob(:,:,:) = 0.0
	t = thresholds_hist
	DO jy = 1, ny
		
		!IF (mod(jy,100) .eq. 0) PRINT *,'      jy = ',jy,' of ',ny
		DO ix = 1, nx
	!DO jy = jmax, jmax
	!	DO ix = imax,imax !  imax -100, imax+100, 20
			
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
			
				! ---- Compute final probability as a weighted sum of
				!      kernel densities exceeding the threshold
			
				ens1d(:) = precip_ens(:,jy,ix)
				CALL quicksort(ens1d, 1, nmembers)
				DO imem = 1, nmembers
					IF (ens1d(imem) .gt. 0.0) THEN
						IF (imem .eq. 1) THEN
							pmean = b0_mean_lowrank + ens1d(imem)*b1_mean_lowrank
							!pstd = b0_std_lowrank + ens1d(imem)*b1_std_lowrank
							pstd = MAX(b0_std_lowrank + ens1d(imem)*b1_std_lowrank,0.01)
						ELSE IF (imem .eq. nmembers) THEN
							pmean = b0_mean_highrank + ens1d(imem)*b1_mean_highrank
							!pstd = b0_std_highrank + ens1d(imem)*b1_std_highrank
							pstd = MAX(b0_std_highrank + ens1d(imem)*b1_std_highrank,0.01)
						ELSE	
							pmean = b0_mean_midrank + ens1d(imem)*b1_mean_midrank
							!pstd = b0_std_midrank + ens1d(imem)*b1_std_midrank
							pstd = MAX(b0_std_midrank + ens1d(imem)*b1_std_midrank, 0.01)
						ENDIF
						DO ithresh = 1, nthresholds
							zscore = (pmean - thresholds(ithresh)) / pstd
							idx = 0
							IF (zscore .lt. -5.0) THEN
								cum = 0.0
							ELSE IF (zscore .ge. -5.0 .and. zscore .le. 5.0) THEN
								idx = INT(zscore*1000)
								cum = cum_table(idx)
							ELSE 
								cum = 1.0
							END IF
							CALL cumnor(zscore, cum, ccum) 
							prob(ithresh, jy,ix) = prob(ithresh, jy,ix) + &
								cum * closest_hist(imem,icat)
						END DO
					END IF
					
					
				END DO			
			END IF

		END DO
	END DO
	PRINT *,'  prob(:,jmax,imax) = ', prob(:,jmax,imax)
	!PRINT *,'  POP prob(jmax,1:nx:10) = ', prob(1,jmax,1:nx:10) 
	
	RETURN
	
	END SUBROUTINE weighted_probability_dressed_f90
