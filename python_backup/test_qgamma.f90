PROGRAM test_qgamma
!  gfortran -Wall -O4 -o test_qgamma.x test_qgamma.f90 qgamma.f90 dgamma.f90 exparg.f90 pgamma.f90 qnorm.f90 pnorm.f90 dnorm.f90 gamma.f90 ipmpar.f90

!cumgam.f90 gamma_inc.f90 error_f.f90 error_fc.f90 exparg.f90 gam1.f90 ipmpar.f90 pgamma.f90 dgamma.f90 
!qgamma.f90 rlog.f90 rexp.f90 dnorm.f90 pnorm.f90 qnorm.f90 gamma.f90 splev.f fpbspl.f

DOUBLE PRECISION :: alpha, beta, gscale, gcum
REAL*8 qmp, qgamma

LOGICAL tootrue, toofalse

tootrue = .false.
toofalse = .false.

alpha = 2.0
beta = 1.0
gscale = 1.0 / beta
PRINT *, gscale
gcum = 0.5
PRINT *,'gcum, alpha, beta, gscale before = ', gcum, alpha, beta, gscale
qmp = qgamma(gcum, alpha, gscale, tootrue, toofalse)
PRINT *,'qmp after = ', qmp

alpha = 2.0   
beta = 10.0    
gscale = 1.0 / beta
gcum = 0.5
PRINT *,'gcum, alpha, beta, gscale before = ', gcum, alpha, beta, gscale
qmp = qgamma(gcum, alpha, gscale, tootrue, toofalse)
PRINT *,'qmp after = ', qmp

STOP
END


