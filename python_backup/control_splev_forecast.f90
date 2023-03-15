SUBROUTINE control_splev_forecast(knots_in, bspline_coef_in, precip_in, nk, &
	nz, ier, cdf_out)
	
! sudo f2py --opt='-O4' --opt='-Wno-tabs' -c -m control_splev_forecast control_splev_forecast.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f
	
INTEGER, INTENT(IN) :: nk, nz
REAL*8, INTENT(IN), DIMENSION(nk) :: knots_in, bspline_coef_in
REAL*8, INTENT(IN), DIMENSION(nz) :: precip_in
INTEGER, INTENT(OUT) :: ier
REAL*8, INTENT(OUT), DIMENSION(nz) :: cdf_out

! f2py intent(in) nk, nz
! f2py intent(in) knots_in, bspline_coef_in, hazard_fn_in
! f2py depend(nk) knots_in, bspline_coef_in
! f2py depend(nz) hazard_fn_in, cdf_out
! f2py intent(out) ier, cdf_out

REAL*4, DIMENSION(nk) :: knots, bspline_coef
REAL*4, DIMENSION(nz) :: precip
REAL*4, DIMENSION(nz) :: cdf
INTEGER nthree

knots = knots_in  ! convert to real*4's as splev needs
bspline_coef = bspline_coef_in
precip = precip_in
nthree = 3

print *,'calling splev'
print *,'  knots = ',knots
print *,'  nk, nthree, nz = ',nk, nthree, nz
CALL splev(knots, nk, bspline_coef, nthree, precip, cdf, nz, ier)
cdf_out = cdf
print *,'back from splev'

RETURN
END SUBROUTINE control_splev_forecast
