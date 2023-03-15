# generate an estimate of pi using the infinite series
# truncated to n terms

def main():

    import math
    n = input("Enter the number of terms in calculation of pi : ")
    sign = -1.0
    sumpi = 0.
    for i in range(n):
        denom = float(i*2)+1.
        sign = sign*-1.0
        sumpi = sumpi + sign*4.0/denom
    print 'estimate, real, difference = ',\
        sumpi,math.pi,abs(sumpi-math.pi)

main()
