SUBROUTINE quantile_mapping_gamma_mixture_v2_f90(&
    spline_info_mswep, indices_to_query_mswep, fzero_mswep, &
    lons_mswep, lats_mswep, spline_info_gefsv12, &
	indices_to_query_gefsv12, fzero_gefsv12, precip_realtime, &
	lons_1d_gefsv12, lats_1d_gefsv12, ncomponents, ny_mswep, &
	nx_mswep, nspline1, nspline2, nspline3, ny_gefsv12, &
	nx_gefsv12, qmapped_precip)

! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m quantile_mapping_gamma_mixture_v2_f90 quantile_mapping_gamma_mixture_v2_f90.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90

INTEGER, INTENT(IN) :: ncomponents, ny_mswep, nx_mswep, ny_gefsv12, &
	nx_gefsv12, nspline1, nspline2, nspline3

REAL*8, INTENT(IN), DIMENSION(ny_mswep, nx_mswep, nspline1, nspline2) :: spline_info_mswep
REAL*8, INTENT(IN), DIMENSION(ny_mswep, nx_mswep, nspline3) :: indices_to_query_mswep
REAL*8, INTENT(IN), DIMENSION(ny_mswep, nx_mswep) :: fzero_mswep
REAL, INTENT(IN), DIMENSION(ny_mswep, nx_mswep) :: lats_mswep, lons_mswep

REAL*8, INTENT(IN), DIMENSION(ny_gefsv12, nx_gefsv12, nspline1, nspline2) :: spline_info_gefsv12
REAL*8, INTENT(IN), DIMENSION(ny_gefsv12, nx_gefsv12, nspline3) :: indices_to_query_gefsv12
REAL*8, INTENT(IN), DIMENSION(ny_gefsv12, nx_gefsv12) :: fzero_gefsv12

REAL, INTENT(IN), DIMENSION(ny_gefsv12) :: lats_1d_gefsv12 
REAL, INTENT(IN), DIMENSION(nx_gefsv12) :: lons_1d_gefsv12

REAL, DIMENSION(ny_mswep, nx_mswep), INTENT(IN) :: precip_realtime

REAL*8, DIMENSION(ny_mswep, nx_mswep), INTENT(OUT) :: qmapped_precip

REAL*8, DIMENSION(10) :: pamounts, anal_quantiles
REAL*8 gefsv12_quantile, conv_criterion, f, frac, &
	pgefs, qdiff, qest, plow, phi, qlow, qhi, &
	estimated_qmapped_precip, estimated_quantile
REAL*8, DIMENSION(ny_mswep, nx_mswep) :: precip_realtime_double	

REAL rlon, rlat

DATA pamounts /0.0001, 0.1, 0.5, 1.0, 3.0,   5.0, 10.0, 25.0, 50.0, 250.0/

! f2py intent(in) ncomponents, ny_mswep, nx_mswep, ny_gefsv12, nx_gefsv12
! f2py intent(in) nspline1, nspline2, nspline3

! f2py intent(in) spline_info_mswep, indices_to_query_mswep
! f2py intent(in) fzero_mswep, lons_mswep, lats_mswep
! f2py depend(ny_mswep, nx_mswep, nspline1, nspline2) spline_info_mswep  
! f2py depend(ny_mswep, nx_mswep, nspline3) indices_to_query_mswep  
! f2py depend(ny_mswep, nx_mswep) fzero_mswep, lats_mswep, lons_mswep

! f2py intent(in) spline_info_gefsv12, indices_to_query_gefsv12
! f2py intent(in) fzero_gefsv12, precip_realtime
! f2py intent(in) lons_1d_gefsv12, lats_1d_gefsv12

! f2py depend(ncomponents, ny_gefsv12, nx_gefsv12, nspline1, nspline2) spline_info_gefsv12
! f2py depend(ncomponents, ny_gefsv12, nx_gefsv12, nspline3) indices_to_query_gefsv12

! f2py depend(ny_gefsv12, nx_gefsv12) fzero_gefsv12
! f2py depend(ny_mswep, nx_mswep) precip_realtime
! f2py depend(ny_gefsv12) lats_1d_gefsv12
! f2py depend(nx_gefsv12) lons_1d_gefsv12

! f2py intent(out) qmapped_precip
! f2py depend(ny_mswep, nx_mswep) qmapped_precip


INTEGER, DIMENSION(ny_mswep, nx_mswep) :: iclose, jclose
LOGICAL readclosest
CHARACTER*80 outfile, infile

conv_criterion = 0.0001
precip_realtime_double = precip_realtime

! ---- determine the GEFSv12 indices of closest grid point to each MSWEP point

readclosest = .false.
IF (readclosest .eqv. .false.) THEN
	PRINT *, 'Determining closest GEFSv12 grid points'
	DO jy = 1, ny_mswep
		IF (jy/10 .eq. 0) THEN
			PRINT *,'Processing jy ',jy,' of ',ny_mswep
		ENDIF
		DO ix = 1, nx_mswep
			rlon = lons_mswep(jy,ix)
			frac = (rlon-lons_1d_gefsv12(1)) / &
				(lons_1d_gefsv12(nx_gefsv12) - lons_1d_gefsv12(1))
			idxx = 1 + INT(REAL(nx_gefsv12)*frac)
			
			IF (idxx .lt. 1 .or. idxx .gt. nx_gefsv12) THEN
				IF (idxx .lt. 1) idxx = 1
				IF (idxx .gt. nx_gefsv12) idxx = nx_gefsv12
			ENDIF
			
			iclose(jy,ix) = idxx
		
			rlat = lats_mswep(jy,ix)
			frac = (rlat-lats_1d_gefsv12(1)) / &
				(lats_1d_gefsv12(ny_gefsv12) - lats_1d_gefsv12(1))
			idxy = 1 + INT(REAL(ny_gefsv12)*frac)
			
			IF (idxy .lt. 1 .or. idxy .gt. ny_gefsv12) THEN
				IF (idxy .lt. 1) idxy = 1
				IF (idxy .gt. ny_gefsv12) idxy = ny_gefsv12
			ENDIF
			
			jclose(jy,ix) = idxy
		END DO
	END DO
	
	outfile = 'ijclosest.dat'
	OPEN (unit=1,file=outfile, status='replace', form='unformatted')
	WRITE (1) iclose
	WRITE (1) jclose
	CLOSE (1)
ELSE
	outfile = 'ijclosest.dat'
	OPEN (unit=1,file=outfile, status='old', form='unformatted')
	READ (1) iclose
	READ (1) jclose
	CLOSE (1)
ENDIF

DO jy = 1, ny_mswep
!DO jy = ny_mswep/2, ny_mswep/2
	
	IF (MOD(jy,10) .eq. 0) PRINT *,'********** jy = ',jy,'    ***********'
	DO ix = 1, nx_mswep
	!	DO ix = nx_mswep/10, 9*nx_mswep/10, 20
		
		!PRINT *,'************************ jy, ix = ',jy, ix, '    ***********'
		
		! ---- estimate the quantile associated with several analyzed  
		!      precipitation amounts 
		
		DO iamt = 1, 10
			CALL get_quantile_mswep(pamounts(iamt), jy, ix, &
				ny_mswep, nx_mswep, nspline1, nspline2, nspline3, &
				spline_info_mswep, indices_to_query_mswep, fzero_mswep, &
				anal_quantiles(iamt))
		END DO
		!PRINT *,'   pamounts: ',  pamounts(:)
		!PRINT *,'   analysis quantiles: ', anal_quantiles(:)
		
		! ---- get the quantile associated with today's GEFSv12 forecast
		!      at the nearest GEFSv12 grid point
		
		ixgefs = iclose(jy,ix)
		jygefs = jclose(jy,ix)
		!print *,'   ixgefs, jygefs = ', ixgefs, jygefs
		pgefs = precip_realtime_double(jy,ix)
		CALL get_quantile(pgefs, jygefs, &
			ixgefs, ny_gefsv12, nx_gefsv12, weights_gefsv12, &
			alpha_gefsv12, beta_gefsv12, fzero_gefsv12, &
			gefsv12_quantile)
		!PRINT *,'   todays forecast = ', pgefs
		!PRINT *,'   gefsv12 quantile = ', gefsv12_quantile
		
		! ---- determine the indices into the anal_quantile array
		!      that bound the gefsv12 quantile	
		
		IF (gefsv12_quantile < anal_quantiles(1)) THEN
			qmapped_precip(jy,ix) = 0.0

		ELSE	
			CALL get_bounding_indices(gefsv12_quantile, 10, &
				anal_quantiles, ilow, ihigh)
			!PRINT *,'   bounding indices ilow, ihigh = ', ilow, ihigh
			!PRINT *,'   gefsv12_quantile, anal low, hi = ', &
			!	gefsv12_quantile, anal_quantiles(ilow), anal_quantiles(ihigh)
			f = (gefsv12_quantile - anal_quantiles(ilow)) / &
				(anal_quantiles(ihigh) - anal_quantiles(ilow))
			!PRINT *,'   f = ', f
			qdiff = 1.0
			qlow = anal_quantiles(ilow) 
			qhi = anal_quantiles(ihigh)
			qest = (1.0 - f)*qlow + f*qhi
			plow = pamounts(ilow)
			phi = pamounts(ihigh)
			!PRINT *,'   qdiff, qest, plow, phi = ', qdiff, qest, plow, phi
			!PRINT *,'   DO WHILE loop!'
			DO WHILE (qdiff .gt. conv_criterion)
				estimated_qmapped_precip = (1.0 - f)*plow + f*phi
				!PRINT *,'   estimated_qmapped_precip = ', estimated_qmapped_precip
				CALL get_quantile(estimated_qmapped_precip, jy, ix, &
					ny_mswep, nx_mswep, weights_mswep, &
					alpha_mswep, beta_mswep, fzero_mswep, &
					estimated_quantile)
				!PRINT *,'   estimated_quantile = ', estimated_quantile
				IF (estimated_quantile .ge. gefsv12_quantile) THEN
					phi = estimated_qmapped_precip
					qhi = estimated_quantile
				ELSE
					plow = estimated_qmapped_precip
					qlow = estimated_quantile
				ENDIF
				qdiff = ABS(estimated_quantile - qest)
				qest = estimated_quantile
				f = (gefsv12_quantile - qlow) / (qhi - qlow)
				!PRINT *,'   qdiff, qest, plow, phi, qlow, qhi, f = ', &
				!	qdiff, qest, plow, phi, qlow, qhi, f
			END DO
			qmapped_precip(jy,ix) = estimated_qmapped_precip
		ENDIF
		!PRINT *,'   qmapped_precip(jy,ix) = ', qmapped_precip(jy,ix)
	END DO
END DO			
		
PRINT *,'Done quantile mapping.'	
RETURN
END SUBROUTINE quantile_mapping_gamma_mixture_v2_f90

! ======================================================================

SUBROUTINE get_quantile_mswep(precip_amount, jy, ix, &
	ny, nx, nspline1, nspline2, nspline3, spline_info_mswep, &
	indices_to_query_mswep, fzero, quantile)
	
REAL*8, INTENT(IN):: precip_amount
INTEGER, INTENT(IN) :: jy, ix, ny, nx, nspline1, nspline2, nspline3 
REAL*8, INTENT(IN), DIMENSION(ny, nx, nspline1, nspline2) :: spline_info_mswep
REAL*8, INTENT(IN), DIMENSION(ny, nx, nspline3) :: indices_to_query_mswep
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: fzero
REAL*8, INTENT(OUT) :: quantile
	
REAL*8 cum, nonexceedance_prob, y, cum

!PRINT *,'   inside get_quantile, precip_amount = ', precip_amount
IF (precip_amount .eq. 0.0) THEN
	quantile = 0.0
ELSE	
	
	! ---- this would be very light precip; use fitted Gamma dist
	
	IF (indices_to_query_mswep(jy,ix,1) .eq. -1) THEN	
		y1 = precip_amount / spline_info_mswep(jy,ix,2,1)
	

	CALL cumgam(y1, alpha(1,jy,ix), cum1, nonexceedance_prob1)
	CALL cumgam(y2, alpha(2,jy,ix), cum2, nonexceedance_prob2)
	CALL cumgam(y3, alpha(3,jy,ix), cum3, nonexceedance_prob3)
	!PRINT *, '   cum1,2,3 = ',cum1,cum2, cum3
	!PRINT *, '   fzero = ', fzero(jy,ix)
	quantile = fzero(jy,ix) + (1. - fzero(jy,ix))* &
		(weights(1,jy,ix)*cum1 + &
	 	weights(2,jy,ix)*cum2 + &
	 	weights(3,jy,ix)*cum3)
	!PRINT *, 'quantile = ', quantile
END IF 

RETURN
END SUBROUTINE get_quantile



! ======================================================================

SUBROUTINE get_quantile(precip_amount, jy, ix, &
	ny, nx, weights, alpha, beta, fzero, quantile)
	
REAL*8, INTENT(IN):: precip_amount
INTEGER, INTENT(IN) :: jy, ix, ny, nx
REAL*8, INTENT(IN), DIMENSION(3, ny, nx) :: weights, alpha, beta
REAL*8, INTENT(IN), DIMENSION(ny, nx) :: fzero
REAL*8, INTENT(OUT) :: quantile
	
REAL*8 cum, nonexceedance_prob1, nonexceedance_prob2, &
	nonexceedance_prob3, y1, y2, y3, cum1, cum2, cum3

!PRINT *,'   inside get_quantile, precip_amount = ', precip_amount
IF (precip_amount .eq. 0.0) THEN
	quantile = 0.0
ELSE		
	y1 = precip_amount / beta(1,jy,ix)
	y2 = precip_amount / beta(2,jy,ix)
	y3 = precip_amount / beta(3,jy,ix)
	!PRINT *,'   y1, y2, y3 = ', y1, y2, y3
	!PRINT *,'   alphas = ', alpha(1,jy,ix), alpha(2,jy,ix), alpha(3,jy,ix)
	!PRINT *,'   betas = ', beta(1,jy,ix), beta(2,jy,ix), beta(3,jy,ix)
	!PRINT *,'   weights = ', weights(1,jy,ix), weights(2,jy,ix), weights(3,jy,ix)

	CALL cumgam(y1, alpha(1,jy,ix), cum1, nonexceedance_prob1)
	CALL cumgam(y2, alpha(2,jy,ix), cum2, nonexceedance_prob2)
	CALL cumgam(y3, alpha(3,jy,ix), cum3, nonexceedance_prob3)
	!PRINT *, '   cum1,2,3 = ',cum1,cum2, cum3
	!PRINT *, '   fzero = ', fzero(jy,ix)
	quantile = fzero(jy,ix) + (1. - fzero(jy,ix))* &
		(weights(1,jy,ix)*cum1 + &
	 	weights(2,jy,ix)*cum2 + &
	 	weights(3,jy,ix)*cum3)
	!PRINT *, 'quantile = ', quantile
END IF 

RETURN
END SUBROUTINE get_quantile


! ======================================================================

SUBROUTINE get_bounding_indices(gefsv12_quantile, nquants, &
	anal_quantiles, ilow, ihigh)
	
REAL*8, INTENT(IN) :: gefsv12_quantile
REAL*8, INTENT(IN), DIMENSION(nquants) :: anal_quantiles
INTEGER, INTENT(OUT) :: ilow, ihigh

LOGICAL foundit

foundit = .false.
ilow = 1
ihigh = 2
DO WHILE (foundit .eqv. .false. .and. ilow .ge. nquants-1)
	IF (gefsv12_quantile .ge. anal_quantiles(ilow) &
	.and. gefsv12_quantile .le. anal_quantiles(ihigh)) THEN
		foundit = .true.
	ELSE
		ilow = ilow+1
		ihigh = ihigh+1
	ENDIF
END DO
!PRINT *,'   ilow, ihigh = ', ilow,ihigh

RETURN
END SUBROUTINE get_bounding_indices

! =======================================================================

tcktuple

    A tuple (t,c,k) containing the vector of knots, the B-spline coefficients, and the degree of the spline.


!************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db1ink]] or one of its
!  derivatives at the point `xval`.
!
!  To evaluate the interpolant itself, set `idx=0`,
!  to evaluate the first partial with respect to `x`, set `idx=1`, and so on.
!
!  [[db1val]] returns 0.0 if (`xval`,`yval`) is out of range. that is, if
!```fortran
!   xval < tx(1) .or. xval > tx(nx+kx)
!```
!  if the knots `tx` were chosen by [[db1ink]], then this is equivalent to:
!```fortran
!   xval < x(1) .or. xval > x(nx)+epsx
!```
!  where
!```fortran
!   epsx = 0.1*(x(nx)-x(nx-1))
!```
!
!  The input quantities `tx`, `nx`, `kx`, and `bcoef` should be
!  unchanged since the last call of [[db1ink]].
!
!### History
!  * Jacob Williams, 10/30/2015 : Created 1D routine.

    pure subroutine db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

    implicit none

    integer(ip),intent(in)               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    integer(ip),intent(in)               :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    real(wp),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction. (same as in last call to [[db1ink]])
    real(wp),dimension(nx),intent(in)    :: bcoef    !! the b-spline coefficients computed by [[db1ink]].
    real(wp),intent(out)                 :: f        !! interpolated value
    integer(ip),intent(out)              :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer(ip),intent(inout)            :: inbvx    !! initialization parameter which must be set
                                                     !! to 1 the first time this routine is called,
                                                     !! and must not be changed by the user.
    real(wp),dimension(3_ip*kx),intent(inout) :: w0  !! work array
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return

    call dbvalu(tx,bcoef,nx,kx,idx,xval,inbvx,w0,iflag,f,extrap)

    end subroutine db1val
		

