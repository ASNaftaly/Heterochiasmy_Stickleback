#calculating enrichment for TSSs covered by peaks from random distributions
#read in permutations file
#need to count total number of unique occurrences at TSSs (so count TSS only once per permutation)
#to run script: python3 Calc_enrichment_TSSs_with_peaks_10kperms.py <permutations file> <chr number> <output>

import sys

#read in permutations file:
#returns dictionary with key == perm number and value == isoform id
def read_perm():
    perm_file = sys.argv[1]
    chr_number = sys.argv[2]
    perm_dict = {}
    with open(perm_file, 'r') as perms:
        for line in perms:
            new_line = line.split()
            perm_number = new_line[0]
            chr_num = new_line[1]
            if chr_num == chr_number:
                isoform_id = new_line[5]
                if perm_number in perm_dict:
                    perm_dict[perm_number].append(isoform_id)
                elif perm_number not in perm_dict:
                    perm_dict.update({perm_number:[isoform_id]})
    return perm_dict

#count isoform ids (no duplicates)
#returns dictionary with key == perm number and value = total isoform counts
def count():
    perms = read_perm()
    perm_totals = {}
    for perm in perms:
        single_perm = perms[perm]
        single_chr_perm_count = list(set(single_perm))
        perm_totals.update({perm:len(single_chr_perm_count)})
    return perm_totals


def write():
    perms = count()
    chr = sys.argv[2]
    output = sys.argv[3]
    with open(output, 'a') as out:
        for p in perms:
            final = "%s\t%s\t%s\n" % (chr, str(p), str(perms[p]))
            out.write(final)


write()
