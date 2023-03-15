! compile for python call using
! f2py -c -m compute_closest_member_dress_beta_f90 compute_closest_member_dress_beta_f90.f90 ran3.f
! f2py --opt='-O4' -c -m compute_closest_member_dress_beta_f90 compute_closest_member_dress_beta_f90.f90 ran3.f

SUBROUTINE compute_closest_member_dress_beta_f90(&
	nmembers_x25, ny, nx, ntwofiftyone, &
	ensmean, precip_ens, precip_anal, &
	sum_fracrank, sum_fracrank_squared, nsamps_fracrank, &
	istat)
	
! following Hamill, T. M., and Scheuerer, M., 2018: Probabilistic precipitation 
! forecast postprocessing using quantile mapping and rank-weighted best-member 
! dressing.  Mon. Wea. Rev., 146, 4079-4098.	
! 
! this procedure tallies for this set of ensembles the closest-member histogram 
! statistics.   These are indexed by the ensemble-mean forecast precipitation
! amount.
!
! Coded by Tom Hamill, July 2021
    
INTEGER*4, INTENT(IN) :: nmembers_x25, ny, nx, ntwofiftyone
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: ensmean
REAL*8, INTENT(IN), DIMENSION(nmembers_x25, ny, nx) :: precip_ens
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: precip_anal
REAL*8, DIMENSION(ntwofiftyone), INTENT(OUT) :: sum_fracrank, sum_fracrank_squared
INTEGER*4, DIMENSION(ntwofiftyone), INTENT(OUT) :: nsamps_fracrank
INTEGER*4, INTENT(OUT) :: istat
    
!f2py intent(in) nmembers_x25, ny, nx, ntwofiftyone
!f2py intent(in) ensmean, precip_ens, precip_anal
!f2py intent(out) sum_fracrank, sum_fracrank_squared, nsamps_fracrank
!f2py intent(out) istat
!f2py depend(ny, nx) ensmean, precip_anal
!f2py depend(ntwofiftyone) sum_fracrank, sum_fracrank_squared, nsamps_fracrank
!f2py depend(nmembers_x25, ny, nx) precip_ens
    
! ---- local variables
       
REAL*8 rclosest, e, diff, a, eclosest, rm, rma, f
INTEGER ibelow, iequal

! --- determine if there is any corrupt ensemble members, and if so, don't 
!     tally stats for this day.  Otherwise proceed.

rm = MINVAL(precip_ens)
rma = MINVAL(precip_anal)

print *,'   max,min precip_ens = ', maxval(precip_ens), minval(precip_ens)
print *,'   max, min anal = ', maxval(precip_anal), minval(precip_anal)
IF (rm .lt. 0.0 .or. rma .lt. -99.999) THEN
    istat = -1
    PRINT *,'   rm, rma = ', rm, rma
    PRINT *,'   identified bad forecast or analysis data, so skip this day.'
ELSE
	istat = 0
    DO ix = 1, nx
        DO jy = 1, ny
		!DO jy = 1, ny
                
            ! --- determine which member is closest to the analyzed and
            !     how many members have values lower than or equal to analyzed
              
            a = precip_anal(jy,ix)
			IF (a .ge. 0.0) THEN
	            iclosest = 1
	            rclosest = 9999.
	            eclosest = 0.0
	            DO imem = 1, nmembers_x25
	                e = precip_ens(imem,jy,ix)
	                diff = ABS(e-a)
	                IF (diff .lt. rclosest .and. e .gt. -99) THEN
	                    rclosest = diff
	                    eclosest = e
	                    iclosest = imem
	                ENDIF
	            END DO

	            ibelow = 0
	            iequal = 0
	            DO imem = 1, nmembers_x25
	                e = precip_ens(imem,jy,ix)
	                IF (imem .eq. iclosest) THEN
	                    continue
	                ELSE
	                    IF (e .lt. eclosest) ibelow = ibelow + 1
	                    IF (e .eq. eclosest) iequal = iequal + 1
	                ENDIF
	            END DO
                
				! --- determine the closest_histogram rank 
				
				IF (iequal .eq. 0) THEN
					iclosest = ibelow + 1			
	            ELSE ! with equals, randomly assign rank
					r = ran3(idum) * iequal
					ir = INT(r)
					IF (ir .gt. iequal) ir = iequal
					iclosest = ibelow + ir + 1
				ENDIF			
				
				! ---- tally up statistics that can be used for estimating dressing
				!      Gaussian statistics.
			    
				idiff_idx = NINT(ensmean(jy,ix)*5.)
				IF (idiff_idx .lt. 1) idiff_idx = 1
				IF (idiff_idx .gt. ntwofiftyone) idiff_idx = ntwofiftyone

				! --- tally a counter of the sum of the fractional rank,
				!     and the fractional rank squared.   This will be used
				!     to later determine beta distribution statistics 
				!     from which closest-member histograms can be inferred.
				
				f = FLOAT(iclosest-1) / FLOAT(nmembers_x25-1)
				sum_fracrank(idiff_idx) = sum_fracrank(idiff_idx) + f
				sum_fracrank_squared(idiff_idx) = sum_fracrank_squared(idiff_idx) + f**2
				nsamps_fracrank(idiff_idx) = nsamps_fracrank(idiff_idx) + 1
				
			ENDIF ! a >= 0
                
        END DO !ix
    END DO ! jy
END IF ! good data

RETURN
END SUBROUTINE compute_closest_member_dress_beta_f90

                    
     
