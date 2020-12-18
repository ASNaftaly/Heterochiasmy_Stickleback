#Determining all classes/families of TEs identified by repeat masker
#Pulls all class and returns only unique classes/families
#to run script: python3 NumTEclasses.py <Repeatmasker file>
#output is to terminal
#Author: Alice Shanfelter, March 2019

import sys

def pull_TEclasses():
    te_file = sys.argv[1]
    te_class = []
    with open(te_file, 'r') as tes:
        for line in tes:
            new_line = line.split()
            #this removes the two header lines, the spacing line between the header and the data; and lines that end with a *
            #lines that end in an * match another entry where that entry has a higher score and at least 80% match with the entry with an *
            #this will make it so all remaining lines have the same length
            if len(new_line) == 0 or new_line[0] == "SW" or new_line[0] == "score" or new_line[(len(new_line)-1)] == "*":
                continue
            else:
                #TE class is element 10 in split line
                te_class.append(new_line[10])
        print(set(te_class))
        print(len(set(te_class)))

pull_TEclasses()
