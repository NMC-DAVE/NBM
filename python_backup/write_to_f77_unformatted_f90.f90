SUBROUTINE write_to_f77_unformatted_f90(outfile, precip_realtime_ens, &
    lons_1d_realtime, lats_1d_realtime,  nmembers, ny, nx, istat)
	
	! compiled as a python-callable routine with 
	! f2py --opt='-O4' -c -m write_to_f77_unformatted_f90 write_to_f77_unformatted_f90.f90
	
	CHARACTER*(*), INTENT(IN) :: outfile
	INTEGER, INTENT(IN) :: nmembers, ny, nx
	REAL, INTENT(IN), DIMENSION (nmembers, ny, nx) :: precip_realtime_ens
	REAL, INTENT(IN), DIMENSION(ny) :: lats_1d_realtime
	REAL, INTENT(IN), DIMENSION(nx) :: lons_1d_realtime
	
	REAL, DIMENSION(nmembers, nx, ny) :: precip_realtime_ens_transposed
	
	! f2py intent(in) outfile, nmembers, ny, nx, precip_realtime_ens
	! f2py intent(in) lats_1d_realtime, lons_1d_realtime
	! f2py depend(nmembers, ny, nx) precip_realtime_ens
	! f2py depend(ny) lats_1d_realtime
	! f2py depend(nx) lons_1d_realtime
	! f2py intent(out) istat
	

	! --- flip from python with a leading j index to fortran with a leading i index
	
	DO j = 1, ny
		DO i = 1, nx
			precip_realtime_ens_transposed(:,i,j) = precip_realtime_ens(:,j,i)
		END DO
	END DO
	
	! --- write to unformatted file.
	
	PRINT *,'writing to ', outfile
	OPEN(unit=1, file=outfile, status='replace', form='unformatted')
	WRITE (1) nx, ny, nmembers 
	WRITE (1) precip_realtime_ens
	WRITE (1) lats_1d_realtime
	WRITE (1) lons_1d_realtime
	CLOSE(1)
	istat = 0
	RETURN
	
	END SUBROUTINE write_to_f77_unformatted_f90