# create_alphanumeric_grade
#    for a number grade 1-5, assign a letter grade A-F

import string  # include string library for the split function.

def main():
    lettergrade =['A','A','B','C','D']
    lettergrade.extend(5*'F')
    #for iter in range(5):
    #    lettergrade.append('F')
    print lettergrade
    print "This program converts a number to letter grade"
    print
    
    # Get the message to encode
    numgrade = input("Please enter the number grade: ")
    numgraded10 = numgrade/10
    print "numgrade/10 = ", numgraded10
    print "The letter grade is:", lettergrade[10-numgraded10]

main()

