# File: chaos2.py
# A simple program illustrating chaotic behavior.

def main():
    print "This program illustrates a chaotic function"
    x1,x2 = input("Enter 2 numbers between 0 and 1: ")
    y = input("How many iterations? ")
    for i in range(y):
        xcombo = 3.9 * xcombo * (1 - xcombo)
        print xcombo

main()
