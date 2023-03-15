# numbers2text.py
#     A program to convert a sequence of ASCII numbers into
#         a string of text.

import string  # include string library for the split function.

def main():
    print "This program converts a sequence of ASCII numbers into"
    print "the string of text that it represents."
    print
    
    # Get the message to encode
    inString = raw_input("Please enter the ASCII-encoded message: ")
    # Loop through each substring and build ASCII message
    character_list = []
    #message = ""
    for numStr in string.split(inString):
        # convert the (sub)string to a number
        asciiNum = eval(numStr)
        # append character to message
        charry = chr(asciiNum)
        character_list.append(charry)
        
    print "The decoded message is:", string.join(character_list,"")

main()
