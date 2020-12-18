#calculate overall enrichment of peaks at TSSs compared to observed values
#first need to total number of peaks at TSSs per autosome per permutation
#read in permutations file
#to run script: python3 Calc_chromosome_total_randompeaks_atTSSs.py <permutations file> <chr number> <output>
#Author: Alice Naftaly, Oct 2020

import sys

#pulling random TSS enrichment distributions from 10k permutations
#input file format:
#permutation number \t chr number \t peak start \t peak end \t tss \t isoform id
#returns dictionary with key = chr number and value = [permutation nuber, peak start, peak end, tss] for each peak/tss in a chromosome
def pull_permutations():
	permutations_file = sys.argv[1]
	chr_num = sys.argv[2]
	permutations_dict = {}
	with open(permutations_file, 'r') as permutations:
		for line in permutations:
			new_line = line.split()
			if new_line[1] == chr_num:
				permutation_number = int(new_line[0])
				peak_start = int(new_line[2])
				peak_end = int(new_line[3])
				tss = int(new_line[4])
				dict_value = [peak_start, peak_end, tss]
				if permutation_number in permutations_dict:
					permutations_dict[permutation_number].append(dict_value)
				elif permutation_number not in permutations_dict:
					permutations_dict.update({permutation_number:[dict_value]})
	return permutations_dict


#totalling peaks at TSSs for each chromosome for each permutation
#returns dictioanry with key == permutation number (0-9999) and value == number of peaks at TSSs for a single permutation
def total_peaks():
	permutations = pull_permutations()
	chr_totals = {}
	for perm in permutations:
		single_perm = permutations[perm]
		total_peaks = len(single_perm)
		chr_totals.update({perm:total_peaks})
	return chr_totals


#write output into a single file
#format:
#chr number \t permutation number \t number of peaks at TSSs
def write():
	counts = total_peaks()
	chr_num = sys.argv[2]
	output = sys.argv[3]
	with open(output, 'a') as out:
		for key in counts:
			single_count = counts[key]
			final = "%s\t%s\t%s\n" % (str(chr_num), str(key), str(single_count))
			out.write(final)

write()
