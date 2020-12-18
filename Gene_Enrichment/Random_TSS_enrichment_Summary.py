#writing summary file for random peak distributions enrichment in TSSs
#summary script should first bin peaks into bins based on TSS location
#will bin one chromosome at a time
#then it should count the total number of peaks in each bin per permutation
#bin number (0-1; in 0.05 increments) total for bin for each permutation
#ex:
#0	394	928	934	.... 10,000 permutations
#1	939	102	22	.... 10,000 permutations
#to run script: python3 Random_TSS_enrichment_Summary.py <random peak TSS enrichment output from Random_ATACpeaks_TSSenrichment_500bp.py> <chr num> <chromosome sizes file> <centromere positions file> <output>
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
				dict_value = [permutation_number, peak_start, peak_end, tss]
				if chr_num in permutations_dict:
					permutations_dict[chr_num].append(dict_value)
				elif chr_num not in permutations_dict:
					permutations_dict.update({chr_num:[dict_value]})
	return permutations_dict


#length of chromosome from read coverage summary files
#format of summary file:
#chr number ("chr" + roman numeral) \t chromosome length
#only autosomes
#returns dictionary: key = chr num, value = chr size
def chr_length():
	chr_size_file = sys.argv[3]
	chr_size_dict = {}
	with open(chr_size_file, 'r') as chr_sizes:
		for line in chr_sizes:
			if line == "Chr\tSize\n":
				continue
			else:
				new_line = line.split()
				chr_num = new_line[0]
				chr_size = new_line[1]
				chr_size_dict.update({chr_num:chr_size})
	return chr_size_dict

#pulls centromere position from file based on chromosome number
#centromere position file format:
#Chr \t Centromere.Pos
#chromosomes are in roman numerals (no "chr" before numeral)
#returns integer of centromere position in bps
def pull_centromere():
	centromere_file = sys.argv[4]
	centromere_dict = {}
	with open(centromere_file, 'r') as centromeres:
		for line in centromeres:
			new_line = line.split()
			if len(new_line) > 0:
				key = new_line[0]
				value = new_line[1]
				centromere_dict.update({key:value})
	return centromere_dict


#create bins for chromosome
#arm 1 represents from the start of the chromosome to the centromere
#arm 2 represents from the centromere to the end of the chromosome
#returns dictionary = key = chr number; value = arm1 bin size, arm 2 bin size
def create_bins():
	chr_size = chr_length()
	centromere_pos = pull_centromere()
	chr_arm_bins_dict = {}
	for chr in chr_size:
		size = chr_size[chr]
		arm_2 = int(size) - int(centromere_pos[chr])
		arm_1 = int(size) - arm_2
		bin_size_arm1 = int(round(arm_1/10,0))
		bin_size_arm2 = int(round(arm_2/10,0))
		#bins variable will have [bin size arm 1, bin size arm 2]
		bin_sizes = [bin_size_arm1,bin_size_arm2]
		chr_arm_bins_dict.update({chr:bin_sizes})
	return chr_arm_bins_dict

#counting number of peaks in bins for each arm
#produces dictionary:
#key = bin number
#value = [chrI bin start, bin end, bin peak count, ...]
def count_permutations_arm1():
	bins = create_bins()
	permutations = pull_permutations()
	all_start = []
	all_end = []
	bin_dict = {}
	for chr in permutations:
		arm_sizes = bins[chr]
		arm1_bin_size = arm_sizes[0]
		single_chr_permutations = permutations[chr]
		for single_peak in single_chr_permutations:
			permutation_count = single_peak[0]
			peak_start = single_peak[1]
			peak_end = single_peak[2]
			all_start.append(peak_start)
			all_end.append(peak_end)
			x = 0
			bin_start = 0
			bin_end = arm1_bin_size
			while x < 10:
				bin_count = 0
				for index, value in enumerate(all_start):
					if int(value) >= bin_start and int(value) <= bin_end and int(all_end[index]) <= bin_end:
						bin_count += 1
					elif int(value) >= bin_start and int(value) <= bin_end and int(all_end[index]) >= bin_end:
						bin_count += 1
					else:
						continue
				if bin_count > 0:
					final = [bin_start, bin_end, permutation_count, bin_count]
					if x in bin_dict:
						bin_dict[x].append(final)
					elif x not in bin_dict:
						bin_dict.update({x:[final]})
				bin_count = 0
				bin_start += arm1_bin_size
				bin_end += arm1_bin_size
				x += 1
			all_start = []
			all_end = []
	return bin_dict

def count_permutations_arm2():
	bins = create_bins()
	centromeres = pull_centromere()
	permutations = pull_permutations()
	chr_num = sys.argv[2]
	all_start = []
	all_end = []
	bin_dict = {}
	for chr in permutations:
		arm_sizes = bins[chr]
		arm2_bin_size = arm_sizes[1]
		centromere_pos = int(centromeres[chr])
		single_chr_permutations = permutations[chr]
		for single_peak in single_chr_permutations:
			permutation_count = single_peak[0]
			peak_start = single_peak[1]
			peak_end = single_peak[2]
			all_start.append(peak_start)
			all_end.append(peak_end)
			x = 10
			bin_start = centromere_pos
			bin_end = arm2_bin_size + centromere_pos
			while x < 20:
				bin_count = 0
				for index, value in enumerate(all_start):
					if int(value) >= bin_start and int(value) <= bin_end and int(all_end[index]) <= bin_end:
						bin_count += 1
					elif int(value) >= bin_start and int(value) <= bin_end and int(all_end[index]) >= bin_end:
						bin_count += 1
					else:
						continue
				if bin_count > 0:
					final = [bin_start, bin_end, permutation_count, bin_count]
					if x in bin_dict:
						bin_dict[x].append(final)
					elif x not in bin_dict:
						bin_dict.update({x:[final]})
				bin_count = 0
				bin_start += arm2_bin_size
				bin_end += arm2_bin_size
				x += 1
			all_start = []
			all_end = []
	return bin_dict


#summarize dictionary
#need to get total counts for a single permutation
def summarize():
	arm1_dict = count_permutations_arm1()
	arm2_dict = count_permutations_arm2()
	chr_num = sys.argv[2]
	combined_dict = {}
	for key in arm1_dict:
		combined_dict.update({key:arm1_dict[key]})
	for key2 in arm2_dict:
		combined_dict.update({key2:arm2_dict[key2]})
	final_dict = {}
	x = 0
	while x < 20:
		if x in combined_dict:
			dict_value = combined_dict[x]
			bin_start = []
			bin_end = []
			perm_number = []
			bin_count = []
			for value in dict_value:
				bin_start.append(value[0])
				bin_end.append(value[1])
				perm_number.append(value[2])
				bin_count.append(value[3])
			set_perm_number = list(set(perm_number))
			y = 0
			while y < 10000:
				single_perm_count = 0
				if y in set_perm_number:
					for index, value in enumerate(perm_number):
						if y == value:
							single_perm_count += bin_count[index]
					if x in final_dict:
						final_dict[x].append(single_perm_count)
					elif x not in final_dict:
						final_dict.update({x:[single_perm_count]})
				elif y not in set_perm_number:
					if x in final_dict:
						final_dict[x].append(0)
					elif x not in final_dict:
						final_dict.update({x:[0]})
				y += 1
		elif x not in combined_dict:
			if x in final_dict:
				final_dict[x].append(0)
			elif x not in final_dict:
				final_dict.update({x:[0]})
		x += 1
	return final_dict


#writing to one output file
#writes all of bin 0, then bin 1, etc.
#no headers will be written
#should be 400 lines per iteration
def write():
	chr_num = sys.argv[2]
	final_dict = summarize()
	bins = ["0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1.0"]
	output = sys.argv[5]
	with open(output,'a') as out:
		x = 0
		while x < 20:
			if x in final_dict:
				dict_value = final_dict[x]
				str_value = [str(i) for i in dict_value]
				joined_value = "\t".join(str_value)
				final = chr_num + "\t" + bins[x] + "\t" + joined_value + "\n"
				out.write(final)
			x += 1


write()
