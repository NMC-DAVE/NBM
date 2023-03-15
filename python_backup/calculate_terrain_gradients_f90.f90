SUBROUTINE calculate_terrain_gradients_f90(earth_radius_meters, &
    terrain_height, lats, ny, nx, terrain_gradx, terrain_grady)
    
! sudo f2py --opt='-O4' -c -m calculate_terrain_gradients_f90 calculate_terrain_gradients_f90.f90
	
INTEGER, INTENT(IN) :: ny, nx
REAL, INTENT(IN) :: earth_radius_meters
REAL, INTENT(IN), DIMENSION (ny, nx) :: terrain_height, lats
REAL*8, INTENT(OUT), DIMENSION (ny, nx) :: terrain_gradx, terrain_grady

! f2py intent(in) earth_radius_meters, terrain_height, lats, ny, nx
! f2py intent(out) terrain_gradx, terrain_grady
! f2py depend(ny,nx) :: terrain_height, lats, terrain_gradx, terrain_grady

rpi = 3.1415926

DO ixa = 1, nxa
	DO jya = 1, nya
    	dy = (2.*rpi*earth_radius_meters) / (360.*2) ! diameter /  1/2 degree
    	dx = dy*cos(lats(jy,1)*rpi/180.)
        
        IF (jy .eq. 1) THEN
            dtdy = (terrain_height(2,ix)-terrain_height(1,ix)) / dy
        ELSE IF (jy .eq. ny) THEN
            dtdy = (terrain_height(ny,ix)-terrain_height(ny-1,ix)) / dy    
        ELSE 
            dtdy = (terrain_height(jy+1,ix)-terrain_height(jy-1,ix)) / (2.*dy) 
        ENDIF
        terrain_grady(jy,ix) = -1. * dtdy
        
        
        IF (ix .eq. 1) THEN
            dtdx = (terrain_height(jy,2)-terrain_height(jy,1)) / dx
        ELSE IF (ix .eq. nx) THEN
            dtdx = (terrain_height(jy,nx)-terrain_height(jy,nx-1)) / dx    
        ELSE 
            dtdx = (terrain_height(jy,ix+1)-terrain_height(jy,ix-1)) / (2.*dx) 
        ENDIF
        
        terrain_gradx(jy,ix) = -1.*dtdx
		
    END DO ! ixa
END DO ! jya
RETURN
END SUBROUTINE calculate_terrain_gradients_f90