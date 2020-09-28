#Calculating expected enrichment of TEs around random peaks
#Will use same logic as TE_observed_counts_ATACpeaks.py
##Will pull ATAC seq peaks and use start and end positions to call peak
#TE can overlap peak in 4 ways (Peak overlaps start of TE, peak overlaps end of TE, peak is in TE, and peak overlaps whole TE)
#Will randomly assign peaks of the same size across chromosomes
#will run separately for each class of TEs to speed things up
#to run script: python3 DNA_TE_broadclass_random_permutations.py <peaks file> <Repeats file> <chromosome sizes file> <output file>
#Author: Alice Naftaly, April 2019, edited Sept 2020

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
#key = feature (DNA, LINES, SINES, LTRs, RC)
#dict values = [chr number, te start pos, te end pos]
def pull_TEs():
    te_file = sys.argv[2]
    te_dict = {}
    with open(te_file, 'r') as te:
        for line in te:
            if len(line) != 124 and len(line) != 131 and len(line) != 1:
                new_line = line.split()
                chr_num = new_line[4]
                feature = new_line[10]
                if feature.startswith("DNA"):
                    final_feature = "DNA"
                    start = new_line[5]
                    end = new_line[6]
                    dict_value = [chr_num, start, end]
                    if final_feature in te_dict:
                        te_dict[final_feature].append(dict_value)
                    elif final_feature not in te_dict:
                        te_dict.update({final_feature:[dict_value]})
    return te_dict

#creating a dictionary of chromosomes and sizes for each sample
#key = chr number (chr roman numeral)
#dictionary value = size of chromosome
def size_dict():
    size_file = sys.argv[3]
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
#key = Chr number
#dictionary value = number of TEs near peaks (no additional space)
def compare_DNA():
    peaks = randomize_peaks()
    TEs = pull_TEs()
    final_te_dict = {}
    for value in TEs["DNA"]:
        chr_num = value[0]
        te_start_final = value[1]
        te_end_final = value[2]
        final = [te_start_final,te_end_final]
        if chr_num in final_te_dict:
            final_te_dict[chr_num].append(final)
        elif chr_num not in final_te_dict:
            final_te_dict.update({chr_num:[final]})
    peak_overlaps_te_start = 0
    peak_overlaps_te = 0
    peak_in_te = 0
    peak_overlaps_te_end = 0
    DNA_TE_dict = {}
    for key in peaks:
        if key in final_te_dict:
            peaks_value = peaks[key]
            te_values = final_te_dict[key]
            chr_num = key
            for val in te_values:
                te_start = int(val[0])
                te_end = int(val[1])
                for v in peaks_value:
                    peak_start = int(v[0])
                    peak_end = int(v[1])
                    #pit
                    if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
                        peak_in_te += 1
                    #pots
                    elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
                        peak_overlaps_te_start += 1
                    #pote
                    elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
                        peak_overlaps_te_end += 1
                    #pot
                    elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
                        peak_overlaps_te += 1
            final_chr_TE_enrichment = [peak_in_te, peak_overlaps_te_start, peak_overlaps_te_end, peak_overlaps_te]
            DNA_TE_dict.update({chr_num:final_chr_TE_enrichment})
            peak_overlaps_te_start = 0
            peak_overlaps_te = 0
            peak_in_te = 0
            peak_overlaps_te_end = 0
    return DNA_TE_dict

#combining final counts and writing them to a file
#final file format:
#header = Chromosome number Peak.in.TE Peak.Overlaps.TE.start Peak.Overlaps.TE.end Peak.overlaps.TE DNA.TE.Total
#last line is the total count for each element class
def combine():
    DNA_counts = compare_DNA()
    chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXX","chrXXI"]
    DNA_output = sys.argv[4]
    DNA_final_total = 0
    with open(DNA_output, 'a') as DNA_out:
        for value in chrs:
            single_DNA = DNA_counts[value]
            DNA_total = single_DNA[0] + single_DNA[1] + single_DNA[2] + single_DNA[3]
            DNA_final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(single_DNA[0]),str(single_DNA[1]),str(single_DNA[2]),str(single_DNA[3]),str(DNA_total))
            DNA_out.write(DNA_final)
combine()
