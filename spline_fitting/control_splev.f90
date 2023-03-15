SUBROUTINE control_splev(knots_in, bspline_coef_in, hazard_fn_in, nk, &
	nz, ier, qmp_out)
	
! f2py --opt='-O4' --opt='-Wno-tabs' -c -m control_splev control_splev.f90 cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f

! then cp the .so file to one that strips out the OS name, e.g., 
! cp ./control_splev.cpython-38-darwin.so ./control_splev.so
	
INTEGER, INTENT(IN) :: nk, nz
REAL*8, INTENT(IN), DIMENSION(nk) :: knots_in, bspline_coef_in
REAL*8, INTENT(IN), DIMENSION(nz) :: hazard_fn_in
INTEGER, INTENT(OUT) :: ier
REAL*8, INTENT(OUT), DIMENSION(nz) :: qmp_out

! f2py intent(in) nk, nz
! f2py intent(in) knots_in, bspline_coef_in, hazard_fn_in
! f2py depend(nk) knots_in, bspline_coef_in
! f2py depend(nz) hazard_fn_in, qmp_out
! f2py intent(out) ier, qmp_out

REAL*4, DIMENSION(nk) :: knots, bspline_coef
REAL*4, DIMENSION(nz) :: hazard_fn
REAL*4, DIMENSION(nz) :: qmp
INTEGER nthree

knots = knots_in  ! convert to real*4's as splev needs
bspline_coef = bspline_coef_in
hazard_fn = hazard_fn_in
nthree = 3

CALL splev(knots, nk, bspline_coef, nthree, &
	hazard_fn, qmp, nz, ier)
qmp_out = qmp

RETURN
END SUBROUTINE control_splev
