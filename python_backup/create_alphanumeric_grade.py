# create_alphanumeric_grade
#    for a number grade 1-5, assign a letter grade A-F

import string  # include string library for the split function.

def main():
    lettergrade =['A','B','C','D','F','F']
    print "This program converts a number to letter grade"
    print
    
    # Get the message to encode
    numgrade = input("Please enter the number grade: ")
    print "The letter grade is:", lettergrade[5-numgrade]

main()
