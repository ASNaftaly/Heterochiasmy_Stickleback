#Making 20 bins per chromosome for random distribution of peaks
#1. need to pull centromere positions, random peaks, & chromosome size
#2. determine size of each bin for each chromosome
#3. count peaks for each bin
#4. write number of peaks per bin to one file
#5. repeat 10,000 times
#to run script: Make20Bins_random_atac.py <centromere positions file> <random peaks file from Random_expectation_peak_distribution.py> <chr size summary file from read coverage files> <permutation output>
#Author: Alice Naftaly, July 2019

import sys

#pulls centromere position from file based on chromosome number
#centromere position file format:
#Chr \t Centromere.Pos
#chromosomes are in roman numerals (no "chr" before numeral)
#returns integer of centromere position in bps
def pull_centromere():
	centromere_file = sys.argv[1]
	centromere_dict = {}
	with open(centromere_file, 'r') as centromeres:
		for line in centromeres:
			if line.startswith("Chr"):
				continue
			else:
				new_line = line.split()
				if len(new_line) > 0:
					key = new_line[0]
					value = new_line[1]
					centromere_dict.update({key:value})
	return centromere_dict

#pulls peaks from random distribution
#peak file format = chr number \t peak start \t peak end
#return dictionary with key = chr and value = [peak start, peak end] for all peaks on that chromosome
def pull_peaks():
	peaks_file = sys.argv[2]
	peaks_dict = {}
	with open(peaks_file, 'r') as peaks:
		for line in peaks:
			new_line = line.split()
			chr_num = new_line[0]
			peak_start = new_line[1]
			peak_end = new_line[2]
			dict_value = [peak_start, peak_end]
			if chr_num in peaks_dict:
				peaks_dict[chr_num].append(dict_value)
			elif chr_num not in peaks_dict:
				peaks_dict.update({chr_num:[dict_value]})
	return peaks_dict

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
			new_line = line.split()
			chr_num = new_line[0]
			chr_size = new_line[1]
			chr_size_dict.update({chr_num:chr_size})
	return chr_size_dict


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
		roman = chr.strip("chr")
		arm_2 = int(size) - int(centromere_pos[roman])
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
def count_peaks_arm1():
	bins = create_bins()
	peaks = pull_peaks()
	all_start = []
	all_end = []
	bin_dict = {}
	for key in bins:
		arm_sizes = bins[key]
		arm1_bin_size = arm_sizes[0]
		single_chr_peaks = peaks[key]
		for val in single_chr_peaks:
			peak_start = val[0]
			peak_end = val[1]
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
			final = [bin_start, bin_end, bin_count]
			if x in bin_dict:
				bin_dict[x].append(final)
			elif x not in bin_dict:
				bin_dict.update({x:[final]})
			bin_start += arm1_bin_size
			bin_end += arm1_bin_size
			x += 1
		all_start = []
		all_end = []
	return bin_dict


def count_peaks_arm2():
	bins = create_bins()
	peaks = pull_peaks()
	centromere_positions = pull_centromere()
	all_start = []
	all_end = []
	bin_dict = {}
	for key in bins:
		arm_sizes = bins[key]
		arm2_bin_size = arm_sizes[1]
		single_chr_peaks = peaks[key]
		roman = key.strip("chr")
		centromere_pos = int(centromere_positions[roman])
		for val in single_chr_peaks:
			peak_start = val[0]
			peak_end = val[1]
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
			final = [bin_start, bin_end, bin_count]
			if x in bin_dict:
				bin_dict[x].append(final)
			elif x not in bin_dict:
				bin_dict.update({x:[final]})
			bin_start += arm2_bin_size
			bin_end += arm2_bin_size
			x += 1
		all_start = []
		all_end = []
	return bin_dict

#writing to one output file
#writes all of bin 0, then bin 1, etc.
#no headers will be written
#should be 400 lines per iteration
def write():
	arm1_dict = count_peaks_arm1()
	arm2_dict = count_peaks_arm2()
	final_dict = {}
	for key in arm1_dict:
		final_dict.update({key:arm1_dict[key]})
	for key2 in arm2_dict:
		final_dict.update({key2:arm2_dict[key2]})
	bins = ["0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1.0"]
	output = sys.argv[4]
	with open(output,'a') as out:
		#header = "Bin.Num\tBin.Start\tBin.End\tNumber.Peaks\n"
		x = 0
		while x < 20:
			if x in final_dict:
				dict_value = final_dict[x]
				for value in dict_value:
					final = "%s\t%s\t%s\t%s\n" % (str(x), str(value[0]), str(value[1]), str(value[2]))
					out.write(final)
				x += 1


write()
