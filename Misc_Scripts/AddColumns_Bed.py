#add columns to bed file
#to run script: python3 AddColumns_Bed.py <input bed file> <strand> <output bed file>

import sys

#read in bed file:
def convert_bed():
    bed_file = sys.argv[1]
    strand = sys.argv[2]
    output = sys.argv[3]
    with open(bed_file, 'r') as bed, open(output, 'a') as out:
        for line in bed:
            new_line = line.split()
            final_line = "%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0],new_line[1],new_line[2], new_line[3], new_line[4], str(strand))
            out.write(final_line)



convert_bed()
