#creating random distribution of peaks to determine how often any given peak should fall in a TSS
#This script will be for the genome (including chr 19, the PAR will be done separately)
#will begin by randomly assigning a value for the total number of peaks identified
#then will pull 500bp in either direction of  TSS, and if random "peak" is in this 1kb region, this will be considered overlapping
#to run script: python3 TSS_expected_enrich_ATACpeaks_500bp.py <chromosome size file> <peak file xls> <Sorted transcripts file> <output file>
#Author: Alice Shanfelter 2018; edited April 2019

import sys
import time
import random

#creating a dictionary of chromosomes and sizes for each sample
#key = chr number (chr roman numeral)
#dictionary value = size of chromosome
def size_dict():
    size_file = sys.argv[1]
    size_dict = {}
    with open(size_file, 'r') as size:
        for line in size:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                chr_size = new_line[1]
                size_dict.update({chr_num:int(chr_size)})
    return size_dict


#pull peak start and end positions
#will randomly identify peaks of the same size across the genome
#returns list of peak sizes per one chromosome
def pull_peaks():
    peak_file = sys.argv[2]
    peak_dict = {}
    with open(peak_file,'r') as peaks:
        for line in peaks:
            if line.startswith("chr\t"):
                continue
            elif line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                start_pos = int(new_line[1])
                end_pos = int(new_line[2])
                if chr_num in peak_dict:
                    peak_dict[chr_num].append(start_pos)
                elif chr_num not in peak_dict:
                    peak_dict.update({chr_num:[start_pos]})
    return peak_dict

#simulates random peaks based on size of chromosome and the number of peaks identified from each chromosome
#dictionary format:
#key = chr number
#value = list of random "peak" midpoints
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
                #print(single_peak_value)
                if chr_num in random_peaks:
                    random_peaks[chr_num].append(random_peak)
                elif chr_num not in random_peaks:
                    random_peaks.update({chr_num:[random_peak]})
                x += 1
        random_peak_list = []
    return random_peaks

#pulls TSSs from Sorted_transcripts.txt file (psl format)
#chr num = 13
#TSS = 15
#returns dictionary in format:
#key = chr number
#value = tss site as an integer
def pull_tss():
    tss_file = sys.argv[3]
    tss_dict = {}
    with open(tss_file, 'r') as tss:
        for line in tss:
            new_line = line.split()
            chr_num = new_line[13]
            tss_start = new_line[15]
            if chr_num in tss_dict:
                tss_dict[chr_num].append(int(tss_start))
            elif chr_num not in tss_dict:
                tss_dict.update({chr_num:[int(tss_start)]})
    return tss_dict

#comparing peaks and TSSs
#returns a dictionary with the format:
#key = chr number
#value = number of TSSs in peaks
def compare_peaks():
    tss = pull_tss()
    peaks = randomize_peaks()
    enriched_peaks = {}
    chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
    for chr in chrs:
        single_chr_tss = tss[chr]
        single_chr_peaks = peaks[chr]
        x = 0
        tss_counts = 0
        while x < len(single_chr_peaks):
            for v in single_chr_tss:
                upstream_tss = int(v) - 500
                downstream_tss = int(v) + 500
                if upstream_tss <= int(single_chr_peaks[x]) <= downstream_tss:
                    tss_counts += 1
            x += 1
        enriched_peaks.update({chr:tss_counts})
        tss_counts = 0
    return enriched_peaks


#writing data to file
def write():  
    final_values = compare_peaks()
    output = sys.argv[4]
    chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
    total_enrich = 0
    with open(output, 'a') as out:
        for chr_num in chrs:
            tss_enrich = final_values[chr_num]
            total_enrich += tss_enrich
            final = "%s\t" % (tss_enrich)
            out.write(final)
        final = "%s\t" % (str(total_enrich))
        out.write(final)
        out.write("\n")

write()
