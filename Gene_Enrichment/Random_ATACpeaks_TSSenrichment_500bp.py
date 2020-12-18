#creating random distribution of peaks to determine how often any given peak should fall in a TSS
#will begin by randomly assigning a value for the total number of peaks identified
#then will pull 500bp in either direction of  TSS, and if random "peak" is in this 1kb region, this will be considered overlapping
#then will write these to a file with:
#Permutation# Chr num Peak start   Peak end TSS Isoform.ID
#to run script: python3 Random_ATACpeaks_TSSenrichment_500bp.py <chromosome size file> <peak file xls> <all isoforms bed file V5 genome <number of permutations> <output file>
#Author: Alice Naftaly, August 2019, edited Sept 2020

import sys
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
                dict_value = [start_pos, end_pos]
                if chr_num in peak_dict:
                    peak_dict[chr_num].append(dict_value)
                elif chr_num not in peak_dict:
                    peak_dict.update({chr_num:[dict_value]})
    return peak_dict

#simulates random peaks based on size of chromosome and the number of peaks identified from each chromosome
#dictionary format:
#key = chr number
#value = [peak start, peak end] for each random peak
def randomize_peaks():
    chr_sizes = size_dict()
    peaks_dict = pull_peaks()
    random_peak_list = []
    peak_size = []
    random_peaks = {}
    for chr_num in chr_sizes:
        single_chr_size = chr_sizes[chr_num]
        single_chr_peaks = peaks_dict[chr_num]
        for value in single_chr_peaks:
            peak_start = int(value[0])
            peak_end = int(value[1])
            peak_size.append(peak_end-peak_start)
        number_single_peaks = len(single_chr_peaks)
        x = 0
        while x < number_single_peaks:
            random_peak = random.randrange(0,single_chr_size)
            if random_peak in random_peak_list:
                continue
            else:
                random_peak_start = random_peak
                random_peak_end = int(random_peak) + peak_size[x]
                random_peak_full = [random_peak_start, random_peak_end]
                random_peak_list.append(int(random_peak))
                if chr_num in random_peaks:
                    random_peaks[chr_num].append(random_peak_full)
                elif chr_num not in random_peaks:
                    random_peaks.update({chr_num:[random_peak_full]})
                x += 1
        random_peak_list = []
    return random_peaks


#pulling TSSs from bed file
#returns list of transcription start sites
def pull_TSS():
    tss_file = sys.argv[3]
    transcript_starts = {}
    with open(tss_file,'r') as tss:
        for line in tss:
            new_line = line.split()
            chrnum = new_line[0]
            isoform_id = new_line[3]
            strand = new_line[5]
            if strand == "+":
                tss = int(new_line[1])
            elif strand == "-":
                tss = int(new_line[2])
            list_value = [tss, isoform_id]
            if chrnum in transcript_starts:
                transcript_starts[chrnum].append(list_value)
            elif chrnum not in transcript_starts:
                transcript_starts.update({chrnum:[list_value]})
    return transcript_starts

#comparing peaks and TSSs
#returns a dictionary with the format:
#key = chr number
#value = number of TSSs in peaks
def compare_peaks():
    tss = pull_TSS()
    peaks = randomize_peaks()
    peak_start = []
    peak_end = []
    enriched_peaks = {}
    chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
    for chr in chrs:
        single_chr_tss = tss[chr]
        single_chr_peaks = peaks[chr]
        for value in single_chr_peaks:
            peak_start.append(value[0])
            peak_end.append(value[1])
        x = 0
        while x < len(single_chr_peaks):
            single_peak_start = peak_start[x]
            single_peak_end = peak_end[x]
            for v in single_chr_tss:
                upstream_tss = int(v[0]) - 500
                downstream_tss = int(v[0]) + 500
                #if peak is in TSS region
                if upstream_tss <= int(single_peak_start) and upstream_tss <= int(single_peak_end) and downstream_tss >= int(single_peak_start) and downstream_tss >= int(single_peak_end):
                    dict_value = [single_peak_start, single_peak_end,v[0],v[1]]
                    if chr in enriched_peaks:
                        enriched_peaks[chr].append(dict_value)
                    elif chr not in enriched_peaks:
                        enriched_peaks.update({chr:[dict_value]})
                #if peak is upstream of TSS, but within 500bp
                elif upstream_tss >= int(peak_start[x]) and upstream_tss < int(peak_end[x]) and downstream_tss >= int(peak_start[x]) and downstream_tss >= int(peak_end[x]):
                    dict_value = [single_peak_start, single_peak_end,v[0],v[1]]
                    if chr in enriched_peaks:
                        enriched_peaks[chr].append(dict_value)
                    elif chr not in enriched_peaks:
                        enriched_peaks.update({chr:[dict_value]})
                #if peak is downstream of TSS, but within 500bp
                elif upstream_tss <= int(peak_start[x]) and upstream_tss <= int(peak_end[x]) and downstream_tss > int(peak_start[x]) and downstream_tss <= int(peak_end[x]):
                    dict_value = [single_peak_start, single_peak_end,v[0],v[1]]
                    if chr in enriched_peaks:
                        enriched_peaks[chr].append(dict_value)
                    elif chr not in enriched_peaks:
                        enriched_peaks.update({chr:[dict_value]})
            x += 1
    return enriched_peaks


#writing data to file
def write():
    final_values = compare_peaks()
    permutation_number = sys.argv[4]
    output = sys.argv[5]
    chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXX","chrXXI"]
    total_enrich = 0
    with open(output, 'a') as out:
        for chr_num in chrs:
            if chr_num not in final_values:
                continue
            else:
                tss_enrich = final_values[chr_num]
                for value in tss_enrich:
                    final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(permutation_number),str(chr_num),str(value[0]),str(value[1]),str(value[2]),str(value[3]))
                    out.write(final)

write()
