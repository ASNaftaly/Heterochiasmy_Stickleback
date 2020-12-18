#Pulling expected counts and calculating expected distribution
#need to exclude chromosome 19
#1. read expeced counts file
#2.pull autosome total counts
#3. total counts across all chromosomes for a permutation and return list of permutations totals
#4. write totals to a file
#to run script: python3 ExpectedDist_atac_TEbroadcounts.py <expected counts file> <output file>
#Author: Alice Naftaly, August 2019

import sys

#read expected counts file
#returns list of total counts per chromosome (should be 210000 values in this list)
def read_expected():
    expected_file = sys.argv[1]
    chr_total = []
    with open(expected_file, 'r') as expected:
        for line in expected:
            if line.startswith("Chr"):
                continue
            else:
                new_line = line.split()
                chr_total.append(new_line[5])
    return chr_total


#calculate total per permutation (should be 20 values per permutation)
#returns list of summed counts per permutation (should be 100000 values in this list)
def calc_total_per_permutation():
    autosomes = read_expected()
    sum_permutations = []
    x = 0
    y = 20
    while x < len(autosomes):
        single_permutation = autosomes[x:y]
        int_single_permutation = [int(i) for i in single_permutation]
        sum_single_permutation = sum(int_single_permutation)
        sum_permutations.append(sum_single_permutation)
        x += 20
        y += 20
    return sum_permutations

#write total permutations to output
#final output format:
#count for a single permutation for each line
def write():
    totals = calc_total_per_permutation()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for value in totals:
            final = "%s\n" % str(value)
            out.write(final)


write()
