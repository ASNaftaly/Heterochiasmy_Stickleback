#total random peaks in TSSs per permutation
#input file format:
#chr num \t permutation number \t number of peaks at TSSs
#to run script: python3 Total_permutations_autosomes_TSS_enrichment.py <output from Calc_chromosome_total_randompeaks_atTSSs.py or Calc_enrichment_TSSs_with_peaks_10kpermutations.py> <output>
#Author: Alice Naftaly, Oct 2020

import sys

#read in permutation file
#returns dictionary with key == perm number and value == number of random peaks at TSSs
def read_permutations():
    perm_file = sys.argv[1]
    perm_dict = {}
    with open(perm_file, 'r') as perms:
        for line in perms:
            new_line = line.split()
            perm_num = new_line[1]
            tss_overlap = int(new_line[2])
            if perm_num in perm_dict:
                perm_dict[perm_num].append(tss_overlap)
            elif perm_num not in perm_dict:
                perm_dict.update({perm_num:[tss_overlap]})
    return perm_dict


#total peaks at TSSs per permutation
#returns dictionary with key == permutation number and value == total random peaks at TSSs across autosomes
def total_perms():
    perms = read_permutations()
    totaled_perms = {}
    for p in perms:
        total_tss_overlap = sum(perms[p])
        totaled_perms.update({p:total_tss_overlap})
    return totaled_perms

#write output
def write():
    perms = total_perms()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for p in perms:
            single_p = perms[p]
            final = "%s\n" % str(single_p)
            out.write(final)


write()
