#calculating expected values for individual TE families to determine which families are enriched (or if all observed are?)
#will run for each broad class of TEs separately, because of running 10k permutations
##Will pull ATAC seq peaks and use start and end positions to call peak
#TE can overlap peak in 4 ways (TE overlaps start of Peak; TE is completely in the peak; the peak is completely in the TE; TE overlaps the end of the peak)
#Will randomly assign peaks of the same size across chromosomes
#to run script: python3 TEfamily_random_expectations_atac.py <peak file.xls> <Repeats file from repeat masker> <feature name (DNA, LINE, SINE, LTR)> <chromosome sizes file> <output>
#Author: Alice Naftaly, July 2019, edited September 2020

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

#pulls TEs from repeat masker file
#format of repeat masker file: SW.repeat.score  Perc.Div    Perc.Del    Perc.Ins    Query.Seq   Query.Start.pos Query.End.pos   Query (left)    +/C     Matching.repeat Repeat.Class/Family     Repeat.Start.pos    Repeat.End.pos  Repeat(left)    ID
#returns dictionary with following format:
#key = feature (DNA, LINES, SINES, LTRs)
#dict values = [chr number, te start pos, te end pos]
def pull_TEs():
	te_file = sys.argv[2]
	feature_name = sys.argv[3]
	te_dict = {}
	with open(te_file, 'r') as te:
		for line in te:
			if len(line) != 124 and len(line) != 131 and len(line) != 1:
				new_line = line.split()
				chr_num = new_line[4]
				feature = new_line[10]
				if feature.startswith(feature_name):
					final_feature = feature_name
					start = new_line[5]
					end = new_line[6]
					te_name = new_line[10]
					dict_value = [chr_num, te_name, start, end]
					if final_feature in te_dict:
						te_dict[final_feature].append(dict_value)
					elif final_feature not in te_dict:
						te_dict.update({final_feature:[dict_value]})
	return te_dict

#creating a dictionary of chromosomes and sizes for each sample
#key = chr number (chr roman numeral)
#dictionary value = size of chromosome
def size_dict():
	size_file = sys.argv[4]
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
#value = peak_start, peak_end
def randomize_peaks():
	chr_sizes = size_dict()
	peaks_dict = pull_peaks()
	random_peak_list = []
	random_peaks = {}
	for chr_num in chr_sizes:
		single_chr_size = chr_sizes[chr_num]
		single_chr_peaks = peaks_dict[chr_num]
		number_single_peaks = len(single_chr_peaks)
		x = 0
		while x < number_single_peaks:
			random_peak = random.randrange(0,single_chr_size)
			if random_peak in random_peak_list:
				continue
			else:
				random_peak_list.append(int(random_peak))
				random_peak_start = random_peak
				random_peak_end = int(random_peak) + single_chr_peaks[x]
				single_peak_value = [random_peak_start, random_peak_end]
				if chr_num in random_peaks:
					random_peaks[chr_num].append(single_peak_value)
				elif chr_num not in random_peaks:
					random_peaks.update({chr_num:[single_peak_value]})
				x += 1
		random_peak_list = []
	return random_peaks


#Compare peaks with each category of TEs
#for each function, the dictionary returned is:
#key = type of te/peak overlap (tops, tip, pit, tope)
#dictionary value = chr_num, te_name, te_start, te_end, peak_start, peak_end
def compare():
	random_peaks = randomize_peaks()
	TEs = pull_TEs()
	final_te_dict = {}
	feature_name = sys.argv[3]
	for value in TEs[feature_name]:
		chr_num = value[0]
		te_name = value[1]
		te_start_final = value[2]
		te_end_final = value[3]
		final = [te_name, te_start_final,te_end_final]
		if chr_num in final_te_dict:
			final_te_dict[chr_num].append(final)
		elif chr_num not in final_te_dict:
			final_te_dict.update({chr_num:[final]})
	TE_random_dict = {}
	for key in random_peaks:
		if key in final_te_dict:
			peaks_value = random_peaks[key]
			te_values = final_te_dict[key]
			chr_num = key
			for val in te_values:
				te_name = str(val[0])
				te_start = int(val[1])
				te_end = int(val[2])
				for v in peaks_value:
					peak_start = int(v[0])
					peak_end = int(v[1])
					#pit
					if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pit", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in TE_random_dict:
							TE_random_dict[chr_num].append(dict_value)
						elif chr_num not in TE_random_dict:
							TE_random_dict.update({chr_num:[dict_value]})
					#pots
					elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pots", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in TE_random_dict:
							TE_random_dict[chr_num].append(dict_value)
						elif chr_num not in TE_random_dict:
							TE_random_dict.update({chr_num:[dict_value]})
					#pote
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pote", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in TE_random_dict:
							TE_random_dict[chr_num].append(dict_value)
						elif chr_num not in TE_random_dict:
							TE_random_dict.update({chr_num:[dict_value]})
					#pot
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pot", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in TE_random_dict:
							TE_random_dict[chr_num].append(dict_value)
						elif chr_num not in TE_random_dict:
							TE_random_dict.update({chr_num:[dict_value]})
	return TE_random_dict

#write output files that show chr number, TE name, te start, te end, peak start peak end for each broad class of TE families
#will have separate output for each broad TE class
def write():
	tes = compare()
	chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXX","chrXXI"]
	output = sys.argv[5]
	feature_name = sys.argv[3]
	with open(output, 'w') as out:
		#total of 7 columns for main part of document
		header = "%sTEs overlapping sigpeaks\nChr\tOverlap.Type\tTE.name\tTE.start\tTE.end\tPeak.Start\tPeak.end\n" % str(feature_name)
		out.write(header)
		for value in chrs:
			if value in tes:
				single_chr = tes[value]
				for val in single_chr:
					final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(val[0]),str(val[1]),str(val[2]),str(val[3]),str(val[4]),str(val[5]))
					out.write(final)
			else:
				continue
write()
