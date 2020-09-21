#pulling TEs based on SW score
#removes any lines with SW score less than 225
#to run script: python3 Filter_TEs_by_SWscore.py <repeat masker .out file> <output file>
#output file is in the same format as the input file
#Author: Alice Naftaly, Sept 2020

import sys

#read repeat masker file
def read_repeat_masker():
    rm_file = sys.argv[1]
    output = sys.argv[2]
    with open(rm_file, 'r') as repeats, open(output, 'a') as out:
        for line in repeats:
            if line.startswith("   SW"):
                header = line
                out.write(header)
            elif len(line) == 148:
                new_line = line.split()
                sw_score = int(new_line[0])
                if sw_score > 225:
                    out.write(line)

read_repeat_masker()
