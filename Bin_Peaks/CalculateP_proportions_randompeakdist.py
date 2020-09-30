#calcualting P value for random peak distribution script
#p value will be for each bin which will provide information on which bins are enriched for peaks or under enriched
#1. will need to read in random peak distribution summary file and observed peak counts for each bin
#2. ask how often is observed count greater than the random values to determien if peaks are enriched in that bin (can I do the opposite to determine which bins are under enriched for peaks?)
#p value calculation = 1 - (#random values less than observed counts/10000)
#4. write p values for each bins (10 total) *will likely need to adjust p value using some correction method due to multiple testing
#to run script: python3 CalculateP_proportions_randompeakdist.py <random peak summary file 10 bins proportions> <observed binned file 10 binns proportions> <output file>
#Author: Alice Naftaly, July 2019, edited November 2019

import sys

#pull random summary
#each line = 10k permutation counts for bins
#returns dictionary with key = bin number and value = 10k permutations
def pull_random_summary():
	permutations_file = sys.argv[1]
	permutations_dict = {}
	with open(permutations_file, 'r') as permutations:
		for line in permutations:
			new_line = line.split()
			bin_key = new_line[0]
			bin_value = new_line[1:len(new_line)]
			permutations_dict.update({bin_key:bin_value})
	return permutations_dict

#pull observed bin counts
#need to condense chromosomes
#returns dictionary with key = bin number and value = observed peak count for all chromosomes
def pull_observed_proportions():
	observed_file = sys.argv[2]
	observed_dict = {}
	with open(observed_file, 'r') as observed:
		for line in observed:
			new_line = line.split()
			if new_line[0] == "Bin.Num":
				continue
			else:
				bin_number = new_line[0]
				bin_proportion = new_line[1]
				if bin_number in observed_dict:
					observed_dict[bin_number].append(bin_proportion)
				elif bin_number not in observed_dict:
					observed_dict.update({bin_number:[bin_proportion]})
	return observed_dict

#comparing observed values with random permutations
#calculates p value for each bin (10 total)
#returns dictionary with key = bin number (1-10) and value = pvalue
def compare():
	permutations = pull_random_summary()
	observed = pull_observed_proportions()
	significance_dict = {}
	for key in observed:
		observed_bin = float(observed[key][0])
		permutation_bin = permutations[key]
		perms_below_observed = 0
		for count in permutation_bin:
			float_count = float(count)
			if float_count <= observed_bin:
				perms_below_observed += 1
		pvalue = round(1 - (perms_below_observed/10000), 6)
		significance_dict.update({key:pvalue})
	return significance_dict

#writing significance values to output
#format of output file:
#Bin number \t Bin Start (relative) \t Bin End (relative) \t Pvalue
def write():
	significance = compare()
	output_file = sys.argv[3]
	bins = ["0","0.05","0.1","0.15","0.2","0.25", "0.3","0.35", "0.4","0.45","0.5"]
	with open(output_file, 'a') as out:
		x = 0
		while x < len(bins)-1:
			single_pvalue = significance[str(bins[x+1])]
			final = "%s\t%s\t%s\t%s\n" % (str(x), str(bins[x]), str(bins[x+1]),str(single_pvalue))
			out.write(final)
			x += 1

write()
