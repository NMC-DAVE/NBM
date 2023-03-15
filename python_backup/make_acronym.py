#  avg_word_length
#    given a sentence, calculate the average word length

import string  # include string library for the split function.

def main():
    lettersum = 0.
    wordcount = 0.
    textin = raw_input("Enter sentence : ")
    for word in string.split(textin):
        numltrs = len(word)
        if word[-1] = "." or
        word[-1] = ";" or
        word[-1] = ":" or
        word[-1] = "!" or
        word[-1] = "," :
            numltrs = numltrs+1
            lettersum = lettersum + numltrs
            wordcount = wordcount + 1.0
    print "The average word length is ',lettersum/wordcount
main()

