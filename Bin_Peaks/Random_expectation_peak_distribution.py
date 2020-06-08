#creating what the random expectation would be if peaks were randomly distributed across the chromosomes
#1. pull observed peak number
#2. randomly call peaks 
#3. sort random peaks into bins (20 per chromosome)
#4. record number of peaks per bin
#5. compare number of peaks per bin with observed values
#to run script: python3 Random_expectation_peak_distribution.py <peak file .xls> <Chromosome_sizes file> <output file>
#Author: Alice Naftaly, July 2019

import sys
import random

#pulls ATAC seq peaks from file (xls) from Macs2
#returns dictionary with following format:
#key = chr number
#dict value = length
def pull_peaks():
	peaks_file = sys.argv[1]
	peaks_dict = {}
	with open(peaks_file, 'r') as peaks:
		for line in peaks:
			if line.startswith("chr"):
				new_line = line.split()
				if new_line[0] == "chr":
					continue
				else:
					chr_num = new_line[0]
					start = int(new_line[1])
					end = int(new_line[2])
					distance = end-start
					if chr_num in peaks_dict:
						peaks_dict[chr_num].append(distance)
					elif chr_num not in peaks_dict:
						peaks_dict.update({chr_num:[distance]})
	return peaks_dict


#creating a dictionary of chromosomes and sizes for each sample
#key = chr number (chr roman numeral)
#dictionary value = size of chromosome
def size_dict():
	size_file = sys.argv[2]
	size_dict = {}
	with open(size_file, 'r') as size:
		for line in size:
			if line.startswith("chr"):
				new_line = line.split()
				chr_num = new_line[0]
				chr_size = new_line[1]
				size_dict.update({chr_num:int(chr_size)})
	return size_dict

#simulates random peaks based on size of chromosome and the number of peaks identified from each chromosome
#dictionary format:
#key = chr number
#value = peak start
def randomize_peaks():
	chr_sizes = size_dict()
	peaks_dict = pull_peaks()
	random_peak_list = []
	random_peaks = {}
	for key in chr_sizes:
		single_chr_size = chr_sizes[key]
		single_chr_peaks = peaks_dict[key]
		number_of_peaks = len(single_chr_peaks)
		x = 0
		while x < number_of_peaks:
			random_peak = random.randrange(0,single_chr_size)
			if random_peak in random_peak_list:
				continue
			else:
				random_peak_list.append(int(random_peak))
				random_peak_start = random_peak
				random_peak_end = int(random_peak) + int(single_chr_peaks[x])
				single_peak_value = [random_peak_start, random_peak_end]
				if key in random_peaks:
					random_peaks[key].append(single_peak_value)
				elif key not in random_peaks:
					random_peaks.update({key:[single_peak_value]})
				x += 1
		random_peak_list = []
	return random_peaks


#write random peaks to file
#output = chr num \t random peak start \t random peak end \n
def write():
	random_peaks = randomize_peaks()
	output_file = sys.argv[3]
	with open(output_file, 'w') as out:
		for key in random_peaks:
			single_value = random_peaks[key]
			for value in single_value:
				start = value[0]
				end = value[1]
				final = "%s\t%s\t%s\n" % (str(key), str(start), str(end))
				out.write(final)


write()