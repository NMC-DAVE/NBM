! compile for python call using
! sudo f2py --opt='-O4' -c -m compute_closest_member_f90 compute_closest_member_f90.f90 ran3.f

SUBROUTINE compute_closest_member_f90(nmembers_x25, ny, nx, nfour, &
    thresh_low, thresh_mod, thresh_high, ensmean, precip_ens, &
    precip_anal, closest_histogram, istat)
	
! following Hamill, T. M., and Scheuerer, M., 2018: Probabilistic precipitation 
! forecast postprocessing using quantile mapping and rank-weighted best-member 
! dressing.  Mon. Wea. Rev., 146, 4079-4098.	
! 
! this procedure tallies for this set of ensembles the closest-member histogram 
! statistics.   These are indexed by the ensemble-mean forecast precipitation
! amount.
!
! Coded by Tom Hamill, July 2021
    
INTEGER, INTENT(IN) :: nmembers_x25, ny, nx, nfour
REAL, INTENT(IN) :: thresh_low, thresh_mod, thresh_high
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: ensmean
REAL*8, INTENT(IN), DIMENSION(nmembers_x25, ny, nx) :: precip_ens
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: precip_anal
INTEGER, INTENT(OUT), DIMENSION(nmembers_x25, nfour) :: closest_histogram
INTEGER, INTENT(OUT) :: istat
    
!f2py intent(in) nmembers_x25, ny, nx, nfour
!f2py intent(in) thresh_low, thresh_mod, thresh_high
!f2py intent(in) ensmean, precip_ens, precip_anal
!f2py intent(out) closest_histogram
!f2py intent(out) istat
!f2py depend(ny, nx) ensmean, precip_anal
!f2py depend(nmembers_x25, ny, nx) precip_ens
!f2py depend(nmembers_x25, nfour) closest_histogram
    
! ---- local variables
       
REAL*8 rclosest, e, diff, a, eclosest, rm, rma
INTEGER ibelow, iequal

! --- determine if there is any corrupt ensemble members, and if so, don't 
!     tally stats for this day.  Otherwise proceed.
    
rm = MINVAL(precip_ens)
rma = MINVAL(precip_anal)
IF (rm .lt. 0.0 .or. rma .lt. 0.0) THEN
    istat = -1
    PRINT *,'rm, rma = ', rm, rma
    PRINT *,'identified bad forecast or analysis data, so skip this day.'
ELSE
	istat = 0
    DO ix = 1, nx
        DO jy = 1, ny
                
            ! --- determine which member is closest to the analyzed and
            !     how many members have values lower than or equal to analyzed
                
            a = precip_anal(jy,ix)
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
                
			! --- increment the histogram counter.
			
            em = ensmean(jy,ix)
			IF (em .lt. thresh_low) THEN
				closest_histogram(iclosest,1) = closest_histogram(iclosest,1) + 1
			ELSE IF (em .ge. thresh_low .and. em .lt. thresh_mod) THEN
				closest_histogram(iclosest,2) = closest_histogram(iclosest,2) + 1
			ELSE IF (em .ge. thresh_mod .and. em .lt. thresh_high) THEN
				closest_histogram(iclosest,3) = closest_histogram(iclosest,3) + 1
			ELSE
				closest_histogram(iclosest,4) = closest_histogram(iclosest,4) + 1
            ENDIF
                
        END DO !ix
    END DO ! jy
END IF ! good data
    
RETURN
END SUBROUTINE compute_closest_member_f90

                    
     
