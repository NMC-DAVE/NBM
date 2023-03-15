# newton's method for calculating square roots

def main():

    import math
    number, iterates = input("Number to find square root, iterates: ")
    guess = number/2
    for i in range(iterates):
        guess = (guess+number/guess)/2.
    print 'guess, actual = ',guess, math.sqrt(number)

main()
