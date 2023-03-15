! compile for python call using
! f2py --opt='-O4' -c -m compute_closest_member_dress_f90 compute_closest_member_dress_f90.f90 ran3.f
!

! this will create an executable callable with python that includes in the name the OS version.
! For example, ./compute_closest_member_dress_f90.cpython-38-darwin.so

! copy this to a generic name, e.g.,
! cp ./compute_closest_member_dress_f90.cpython-38-darwin.so ./compute_closest_member_dress_f90.so

SUBROUTINE compute_closest_member_dress_f90(&
	nmembers_x25, ny, nx, nseven, ntwofiftyone, &
    thresh_low, thresh_lowmod, thresh_mod, &
	thresh_modhigh, thresh_high, thresh_superhigh, &
	ensmean, precip_ens, &
    precip_anal, closest_histogram, &
	sumxi_low, sumxi2_low, nsamps_low, &
	sumxi_mid, sumxi2_mid, nsamps_mid, &
	sumxi_high, sumxi2_high, nsamps_high, &
	istat)
	
! following Hamill, T. M., and Scheuerer, M., 2018: Probabilistic precipitation 
! forecast postprocessing using quantile mapping and rank-weighted best-member 
! dressing.  Mon. Wea. Rev., 146, 4079-4098.	
! 
! this procedure tallies for this set of ensembles the closest-member histogram 
! statistics and dressing statistics.   These are indexed by the ensemble-mean 
! forecast precipitation amount.
!
! Coded by Tom Hamill, July 2021
    
INTEGER, INTENT(IN) :: nmembers_x25, ny, nx, nseven, ntwofiftyone
REAL, INTENT(IN) :: thresh_low, thresh_mod, thresh_high 
REAL, INTENT(IN) :: thresh_lowmod, thresh_modhigh, thresh_superhigh
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: ensmean
REAL*8, INTENT(IN), DIMENSION(nmembers_x25, ny, nx) :: precip_ens
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: precip_anal
INTEGER, INTENT(OUT), DIMENSION(nmembers_x25, nseven) :: closest_histogram
REAL*8, DIMENSION(ntwofiftyone), INTENT(OUT) :: sumxi_low, sumxi2_low
REAL*8, DIMENSION(ntwofiftyone), INTENT(OUT) :: sumxi_mid, sumxi2_mid
REAL*8, DIMENSION(ntwofiftyone), INTENT(OUT) :: sumxi_high, sumxi2_high
INTEGER, DIMENSION(ntwofiftyone), INTENT(OUT) :: nsamps_low, nsamps_mid, nsamps_high
INTEGER, INTENT(OUT) :: istat
    
!f2py intent(in) nmembers_x25, ny, nx, nseven, ntwofiftyone
!f2py intent(in) thresh_low, thresh_mod, thresh_high
!f2py intent(in) thresh_lowmod, thresh_modhigh, thresh_superhigh
!f2py intent(in) ensmean, precip_ens, precip_anal
!f2py intent(out) closest_histogram 
!f2py intent(out) sumxi_low, sumxi2_low, nsamps_low
!f2py intent(out) sumxi_mid, sumxi2_mid, nsamps_mid
!f2py intent(out) sumxi_high, sumxi2_high, nsamps_high
!f2py intent(out) istat
!f2py depend(ny, nx) ensmean, precip_anal
!f2py depend(ntwofiftyone) sumxi_low, sumxi2_low, nsamps_low
!f2py depend(ntwofiftyone) sumxi_mid, sumxi2_mid, nsamps_mid
!f2py depend(ntwofiftyone) sumxi_high, sumxi2_high, nsamps_high
!f2py depend(nmembers_x25, ny, nx) precip_ens
!f2py depend(nmembers_x25, nseven) closest_histogram
    
! ---- local variables
       
REAL*8 rclosest, e, diff, a, eclosest, rm, rma
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
			    
				idiff_idx = NINT(eclosest*10.)
				!print *,'idiff_idx, eclosest, iclosest, rclosest = ',idiff_idx,eclosest, iclosest, rclosest
				!print *,'ntwofiftyone = ', ntwofiftyone
				IF (idiff_idx .ge. 1 .and. idiff_idx .le. ntwofiftyone) THEN
					IF (iclosest .eq. 1) THEN
						sumxi_low(idiff_idx) = sumxi_low(idiff_idx) + a
						sumxi2_low(idiff_idx) = sumxi2_low(idiff_idx) + a**2
						nsamps_low(idiff_idx) = nsamps_low(idiff_idx) + 1
					ELSE IF (iclosest .eq. nmembers_x25) THEN
						sumxi_high(idiff_idx) = sumxi_high(idiff_idx) + a
						sumxi2_high(idiff_idx) = sumxi2_high(idiff_idx) + a**2
						nsamps_high(idiff_idx) = nsamps_high(idiff_idx) + 1	
					ELSE
						sumxi_mid(idiff_idx) = sumxi_mid(idiff_idx) + a
						sumxi2_mid(idiff_idx) = sumxi2_mid(idiff_idx) + a**2
						nsamps_mid(idiff_idx) = nsamps_mid(idiff_idx) + 1	
					ENDIF
				END IF
				
				! --- increment the histogram counter.
			
	            em = ensmean(jy,ix)
				IF (em .lt. thresh_low) THEN
					closest_histogram(iclosest,1) = closest_histogram(iclosest,1) + 1
				ELSE IF (em .ge. thresh_low .and. em .lt. thresh_lowmod) THEN
					closest_histogram(iclosest,2) = closest_histogram(iclosest,2) + 1
				ELSE IF (em .ge. thresh_lowmod .and. em .lt. thresh_mod) THEN
					closest_histogram(iclosest,3) = closest_histogram(iclosest,3) + 1
				ELSE IF (em .ge. thresh_mod .and. em .lt. thresh_modhigh) THEN
					closest_histogram(iclosest,4) = closest_histogram(iclosest,4) + 1
				ELSE IF (em .ge. thresh_modhigh .and. em .lt. thresh_high) THEN
					closest_histogram(iclosest,5) = closest_histogram(iclosest,5) + 1
				ELSE IF (em .ge. thresh_high .and. em .lt. thresh_superhigh) THEN
					closest_histogram(iclosest,6) = closest_histogram(iclosest,6) + 1
				ELSE
					closest_histogram(iclosest,7) = closest_histogram(iclosest,7) + 1
	            ENDIF
			ENDIF ! a >= 0
                
        END DO !ix
    END DO ! jy
END IF ! good data

RETURN
END SUBROUTINE compute_closest_member_dress_f90

                    
     
