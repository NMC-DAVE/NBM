SUBROUTINE calculate_cov(x1, x2, n, cov) 
INTEGER, INTENT(IN) :: n 
REAL*8, INTENT(IN), DIMENSION(n) :: x1, x2
REAL*8, INTENT(OUT) :: cov
REAL*8 denom

denom = n-1
cov = SUM(x1*x2) / denom
RETURN 
END SUBROUTINE calculate_cov