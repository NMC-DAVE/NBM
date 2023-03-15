# generate the nth iterate of fibonacci sequence

def main():

    fib1, fib2 = 1., 1.
    n = input("Enter number in Fibonacci sequence (>2) : ")
    for i in range(2,n):
        fib1, fib2 = fib2, fib1+fib2
    print 'fibonacci for n = ',fib2

main()
