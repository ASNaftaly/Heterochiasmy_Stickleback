#writing summary file for random peak distributions
#summary script should combine each x value (20 per x per permutation) for each permutation into a total value for that x bin for that permutation
#then need to combine total for each bin across all permutations to give the following format:
#bin number (0-1; in 0.05 increments) total for bin for each permutation
#ex:
#0	394	928	934	.... 10,000 permutations
#1	939	102	22	.... 10,000 permutations
#to run script: python3 Random_peak_dist_summary.py <random peak dist output from Make20Bins_random_atac.py> <output>
#Author: Alice Naftaly, July 2019

import sys

#pulling random distributions from 10k permutations
#input file format:
#Bin number \t bin start \t bin end \t bin peak count
#returns dictionary with key = bin number (0-19) and value = [bin count for each chromosome 20]
#order of value = 20 chromosomes from first permutation, 20 chromosomes from second permutation, ...
def pull_permutations():
	permutations_file = sys.argv[1]
	permutations_dict = {}
	with open(permutations_file, 'r') as permutations:
		for line in permutations:
			new_line = line.split()
			bin_number = new_line[0]
			bin_count = new_line[3]
			if bin_number in permutations_dict:
				permutations_dict[bin_number].append(bin_count)
			elif bin_number not in permutations_dict:
				permutations_dict.update({bin_number:[bin_count]})
	return permutations_dict

#combining 20 autosomes for each permutation for each bin size
#returns dictionary with key = bin number as in pull_permutations() and value = total count for each bin for each permutation
def combine_chrs():
	permutations = pull_permutations()
	combined_permutations = {}
	for key in permutations:
		single_permutation = permutations[key]
		chr_count_start = 0
		chr_count_end = 20
		while chr_count_start < 200000:
			all_chrs = single_permutation[chr_count_start:chr_count_end]
			int_all_chrs = [int(i) for i in all_chrs]
			sum_all_chrs = sum(int_all_chrs)
			if key in combined_permutations:
				combined_permutations[key].append(sum_all_chrs)
			elif key not in combined_permutations:
				combined_permutations.update({key:[sum_all_chrs]})
			chr_count_start += 20
			chr_count_end += 20
	return combined_permutations

#write permutations to summary output file
#format:
# bin number (0-1) \t 10k permutation bin totals across all chromosomes
def write():
	permutations = combine_chrs()
	bins = ["0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.0"]
	output_file = sys.argv[2]
	with open(output_file, 'a') as out:
		for index, value in enumerate(bins):
			single_bin = permutations[str(index)]
			str_single_bin = [str(i) for i in single_bin]
			tab_permutations = "\t".join(str_single_bin)
			final = str(value) + "\t" + tab_permutations + "\n"
			out.write(final)

write()
