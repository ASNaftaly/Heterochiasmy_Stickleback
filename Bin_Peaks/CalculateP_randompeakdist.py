#calcualting P value for random peak distribution script
#p value will be for each bin which will provide information on which bins are enriched for peaks or under enriched
#1. will need to read in random peak distribution summary file and observed peak counts for each bin
#2. collapse bins to just be between 0-0.5 by combining as follows: [0,19], [1,18], [2,17], [3,16], [4,15], [5,14], [6,13], [7,12], [8,11], [9,10] = 10 bins total
#3. ask how often is observed count greater than the random values to determien if peaks are enriched in that bin (can I do the opposite to determine which bins are under enriched for peaks?)
#p value calculation = 1 - (#random values less than observed counts/10000) 
#4. write p values for each bins (10 total) *will likely need to adjust p value using some correction method due to multiple testing
#to run script: python3 CalculateP_randompeakdist.py <random peak summary file> <observed binned file> <output file> 
#Author: Alice Naftaly, July 2019S

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


#collapse bins from 20 bins to 10 bins for random peaks
#this could be written more succinctly, but this works
#returns dictionary with key = bin number (1-10) and value = collapses permutations counts for 2 bins that represent the same region of a chromosomal arm
def collapse_permutations_random():
	permutations = pull_random_summary()
	collapsed_perm_dict = {}
	bin_1 = permutations["0.05"]
	bin_2 = permutations["0.10"]
	bin_3 = permutations["0.15"]
	bin_4 = permutations["0.20"]
	bin_5 = permutations["0.25"]
	bin_6 = permutations["0.30"]
	bin_7 = permutations["0.35"]
	bin_8 = permutations["0.40"]
	bin_9 = permutations["0.45"]
	bin_10 = permutations["0.50"]
	bin_11 = permutations["0.55"]
	bin_12 = permutations["0.60"]
	bin_13 = permutations["0.65"]
	bin_14 = permutations["0.70"]
	bin_15 = permutations["0.75"]
	bin_16 = permutations["0.80"]
	bin_17 = permutations["0.85"]
	bin_18 = permutations["0.90"]
	bin_19 = permutations["0.95"]
	bin_20 = permutations["1.0"]
	combined_bins_1_20 = [int(x) + int(y) for x, y in zip(bin_1, bin_20)]
	combined_bins_2_19 = [int(x) + int(y) for x, y in zip(bin_2, bin_19)]
	combined_bins_3_18 = [int(x) + int(y) for x, y in zip(bin_3, bin_18)]
	combined_bins_4_17 = [int(x) + int(y) for x, y in zip(bin_4, bin_17)]
	combined_bins_5_16 = [int(x) + int(y) for x, y in zip(bin_5, bin_16)]
	combined_bins_6_15 = [int(x) + int(y) for x, y in zip(bin_6, bin_15)]
	combined_bins_7_14 = [int(x) + int(y) for x, y in zip(bin_7, bin_14)]
	combined_bins_8_13 = [int(x) + int(y) for x, y in zip(bin_8, bin_13)]
	combined_bins_9_12 = [int(x) + int(y) for x, y in zip(bin_9, bin_12)]
	combined_bins_10_11 = [int(x) + int(y) for x, y in zip(bin_10, bin_11)]
	collapsed_perm_dict.update({"1":combined_bins_1_20})
	collapsed_perm_dict.update({"2":combined_bins_2_19})
	collapsed_perm_dict.update({"3":combined_bins_3_18})
	collapsed_perm_dict.update({"4":combined_bins_4_17})
	collapsed_perm_dict.update({"5":combined_bins_5_16})
	collapsed_perm_dict.update({"6":combined_bins_6_15})
	collapsed_perm_dict.update({"7":combined_bins_7_14})
	collapsed_perm_dict.update({"8":combined_bins_8_13})
	collapsed_perm_dict.update({"9":combined_bins_9_12})
	collapsed_perm_dict.update({"10":combined_bins_10_11})
	return collapsed_perm_dict
	
#pull observed bin counts
#need to condense chromosomes
#returns dictionary with key = bin number and value = observed peak count for all chromosomes
def pull_observed_counts():
	observed_file = sys.argv[2]
	observed_dict_uncollapsed = {}
	observed_dict_collapse = {}
	with open(observed_file, 'r') as observed:
		for line in observed:
			new_line = line.split()
			if new_line[0] == "Bin.Num":
				continue
			else:
				bin_number = new_line[0]
				bin_count = new_line[3]
				if bin_number in observed_dict_uncollapsed:
					observed_dict_uncollapsed[bin_number].append(bin_count)
				elif bin_number not in observed_dict_uncollapsed:
					observed_dict_uncollapsed.update({bin_number:[bin_count]})
	for key in observed_dict_uncollapsed:
		single_key = observed_dict_uncollapsed[key]
		int_single_key = [int(i) for i in single_key]
		sum_single_key = sum(int_single_key)
		observed_dict_collapse.update({key:sum_single_key})
	return observed_dict_collapse


#collapse bins from 20 bins to 10 bins for observed peaks
#this could be written more succinctly, but this works
#returns dictionary with key = bin number (1-10) and value = collapses permutations counts for 2 bins that represent the same region of a chromosomal arm
def collapse_permutations_observed():
	permutations = pull_observed_counts()
	collapsed_perm_dict = {}
	combined_bins_1_20 = int(permutations["0.05"]) + int(permutations["1.0"])
	combined_bins_2_19 = int(permutations["0.1"]) + int(permutations["0.95"])
	combined_bins_3_18 = int(permutations["0.15"]) + int(permutations["0.9"])
	combined_bins_4_17 = int(permutations["0.2"]) + int(permutations["0.85"])
	combined_bins_5_16 = int(permutations["0.25"]) + int(permutations["0.8"])
	combined_bins_6_15 = int(permutations["0.3"]) + int(permutations["0.75"])
	combined_bins_7_14 = int(permutations["0.35"]) + int(permutations["0.7"])
	combined_bins_8_13 = int(permutations["0.4"]) + int(permutations["0.65"])
	combined_bins_9_12 = int(permutations["0.45"]) + int(permutations["0.6"])
	combined_bins_10_11 = int(permutations["0.5"]) + int(permutations["0.55"])
	collapsed_perm_dict.update({"1":combined_bins_1_20})
	collapsed_perm_dict.update({"2":combined_bins_2_19})
	collapsed_perm_dict.update({"3":combined_bins_3_18})
	collapsed_perm_dict.update({"4":combined_bins_4_17})
	collapsed_perm_dict.update({"5":combined_bins_5_16})
	collapsed_perm_dict.update({"6":combined_bins_6_15})
	collapsed_perm_dict.update({"7":combined_bins_7_14})
	collapsed_perm_dict.update({"8":combined_bins_8_13})
	collapsed_perm_dict.update({"9":combined_bins_9_12})
	collapsed_perm_dict.update({"10":combined_bins_10_11})
	return collapsed_perm_dict

#comparing observed values with random permutations
#calculates p value for each bin (10 total)
#returns dictionary with key = bin number (1-10) and value = pvalue 
def compare():
	permutations = collapse_permutations_random()
	observed = collapse_permutations_observed()
	significance_dict = {}
	for key in observed:
		observed_bin = observed[key]
		permutation_bin = permutations[key]
		perms_below_observed = 0
		for count in permutation_bin:
			int_count = int(count)
			if count <= observed_bin:
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
		while x < 10:
			key = x + 1
			single_pvalue = significance[str(key)]
			final = "%s\t%s\t%s\t%s\n" % (str(key), str(bins[x]), str(bins[x+1]),str(single_pvalue))
			out.write(final)
			x += 1

write()					
