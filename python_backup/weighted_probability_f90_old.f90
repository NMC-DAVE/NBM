SUBROUTINE weighted_probability_f90(precip_mean, precip_ens, pthresh, &
	thresholds, closest_hist_cdf, ny, nx, nmembers, nchcats, prob)

	! f2py --opt='-O4' -c -m weighted_probability_f90 weighted_probability_f90.f90 
	
	INTEGER, INTENT(IN) :: ny, nx, nmembers, nchcats
	REAL, INTENT(IN), DIMENSION(ny, nx) :: precip_mean
	REAL, INTENT(IN), DIMENSION(nmembers, ny, nx) :: precip_ens
	REAL, INTENT(IN), DIMENSION(nchcats-1) :: thresholds
	REAL, INTENT(IN), DIMENSION(nmembers, nchcats) :: closest_hist_cdf
	REAL, INTENT(OUT), DIMENSION(ny, nx) :: prob

	!f2py intent(in) ny, nx, nmembers, nchcats
	!f2py intent(in) precip_ens, pthresh
	!f2py intent(in) thresholds, closest_hist_cdf
	!f2py intent(out) prob
	!f2py depend(ny,nx) precip_mean, prob
	!f2py depend(nmembers, ny, nx) precip_ens
	!f2py depend(nchcats-1) thresholds
	!f2py depend(nmembers, nchcats) closest_hist_cdf

	! --- local variables
	
	REAL*8 ens1d(nmembers), ens1d_sorted(nmembers)
	INTEGER bines(nmembers)
	REAL t(nchcats-1)
	
	! ---- loop over grid points
	
	t = thresholds
	DO jy = 1, ny
		DO ix = 1, nx
			
			! ---- determine the category for the closest-member
			!      histogram CDF to apply for this point based
			!      on quantile-mapped ensemble-mean precipitation
			
			p = precip_mean(jy,ix)
			IF (p .le. t(1)) THEN
				icat = 1
			ELSE IF (p .gt. t(1) .and. p .lt. t(2)) THEN
				icat = 2
			ELSE IF (p .gt. t(2) .and. p .lt. t(3)) THEN
				icat = 3
			ELSE
				icat = 4
			ENDIF
			
			! ---- count up the number of ensemble members less than
			!      this threshold.
			
			ens1d(:) = precip_ens(:,jy,ix)
			CALL quicksort(ens1d,1,nmembers)
			WHERE (ens1d > pthresh)
				bines = 1
			ELSEWHERE
				bines = 0
			ENDWHERE
			
			! ---- set the probability as the sum of the weights
			!      up to this number of members with precip below
			!      the threshold.   
			
			!isb = SUM(bines)
			isb = nmembers - SUM(bines)
			
			subroutine quicksort(a, first, last)
			
			
			
			IF (isb .eq. 0) THEN
				prob(jy,ix) = 0.0
			ELSE 
				prob(jy,ix) = closest_hist_cdf(isb,icat)
			ENDIF

		END DO
	END DO

	RETURN
	
	END SUBROUTINE weighted_probability_f90


