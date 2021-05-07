#creating 20 bins across chromosome with 5 bins of each size on each size of the centromere for ATAC-seq (read coverage and peaks)
#this script handles peaks with average recombination rates per peaks
#file format: Chr \t Peak.start \t Peak.end \t average recombination rate
#will first pull chromosome and centromere position and figure out size of bins for each arm
#then will pull peaks and rates in each bin
#to run script: python3 Make20Bins_ATAC_peaks_rates.py <centromere positions file with roman numerals> <chr number as roman numeral> <peaks file> <chr num as integer> <length file file> <output file>
#Author: Alice Naftaly, May 2021

import sys
import time

#pulls centromere position from file based on chromosome number
#returns integer of centromere position in Mb
def pull_centromere():
    centromere_file = sys.argv[1]
    chr_num = sys.argv[2]
    with open(centromere_file, 'r') as centromeres:
        for line in centromeres:
            if line.startswith(chr_num + "\t"):
                centromere_pos = int(line.split()[1])
    return centromere_pos

#pulls peaks from peaks files (either significant or all peaks)
#returns list with [peak 1 start, peak 1 end, peak 2 start, peak 2 end, ...]
def pull_peaks():
    peaks_file = sys.argv[3]
    chr_num_roman = sys.argv[2]
    peaks_dict = {}
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            chr_num_full = chr_num_roman + "\t"
            if line.startswith(chr_num_full):
                new_line = line.split()
                peak_start = int(new_line[1])
                peak_end = int(new_line[2])
                rate = float(new_line[3])
                dict_value = [peak_start, peak_end, rate]
                peaks_dict.update({peak_start:dict_value})
    return peaks_dict

#pulls length of chromosome from readcov file
#returns integer of chromosome length
def pull_chr_length():
    chr_length_file = sys.argv[5]
    chr_num_roman = sys.argv[2]
    chr_length = {}
    with open(chr_length_file, 'r') as readcov_file:
        for line in readcov_file:
            final_line = line.split()
            chr_length.update({final_line[0]:final_line[1]})
        final_chr_length = chr_length[chr_num_roman]
    return final_chr_length

#create bins for chromosome
#arm 1 represents from the start of the chromosome to the centromere
#arm 2 represents from the centromere to the end of the chromosome
def create_bins():
    chr_length = int(pull_chr_length())
    centromere_pos = int(pull_centromere())
    arm_2 = chr_length - centromere_pos
    arm_1 = chr_length - arm_2
    bin_size_arm1 = int(round(arm_1/10,0))
    bin_size_arm2 = int(round(arm_2/10,0))
    #bins variable will have [bin size arm 1, bin size arm 2]
    bin_sizes = [bin_size_arm1,bin_size_arm2]
    return bin_sizes

#counting number of peaks in bins for each arm

def record_peaks_arm1():
    bins = create_bins()
    arm1_bin_size = bins[0]
    peaks = pull_peaks()
    bin_dict = {}
    x = 0
    bin_start = 0
    bin_end = arm1_bin_size
    bin_list = []
    while x < 10:
        for key in peaks:
            single_key = peaks[key]
            peak_start = single_key[0]
            peak_end = single_key[1]
            peak_rate = single_key[2]
            if peak_start >= bin_start and peak_start <= bin_end and peak_end <= bin_end:
                bin_list.append(single_key)
            elif peak_start >= bin_start and peak_start <= bin_end and peak_end >= bin_end:
                bin_list.append(single_key)
            else:
                continue
        bin_dict.update({x:bin_list})
        bin_start += arm1_bin_size
        bin_end += arm1_bin_size
        bin_list = []
        x += 1
    return bin_dict

def record_peaks_arm2():
    bins = create_bins()
    arm2_bin_size = bins[1]
    peaks = pull_peaks()
    centromere_pos = pull_centromere()
    bin_dict = {}
    x = 10
    bin_start = centromere_pos
    bin_end = arm2_bin_size + centromere_pos
    bin_list = []
    while x < 20:
        for key in peaks:
            single_key = peaks[key]
            peak_start = single_key[0]
            peak_end = single_key[1]
            peak_rate = single_key[2]
            if peak_start >= bin_start and peak_start <= bin_end and peak_end <= bin_end:
                bin_list.append(single_key)
            elif peak_start >= bin_start and peak_start <= bin_end and peak_end >= bin_end:
                bin_list.append(single_key)
            else:
                continue
        bin_dict.update({x:bin_list})
        bin_start += arm2_bin_size
        bin_end += arm2_bin_size
        bin_list = []
        x += 1
    return bin_dict

#writing to one output file
def write():
    arm1_dict = record_peaks_arm1()
    arm2_dict = record_peaks_arm2()
    chr_int = sys.argv[4]
    bins = ["0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1.0"]
    output = sys.argv[6]
    with open(output,'a') as out:
        header = "Chr.Num\tBin.Num\tPeak.Start\tPeak.End\tAve.Peak.Recomb.Rate\n"
        out.write(header)
        x = 0
        while x < 20:
            if x in arm1_dict:
                dict_value = arm1_dict[x]
                if len(dict_value) >= 1:
                    for val in dict_value:
                        final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr_int),str(bins[x]), str(val[0]),str(val[1]),str(val[2]))
                        out.write(final)
                x += 1
            elif x in arm2_dict:
                d_value = arm2_dict[x]
                if len(d_value) >= 1:
                    for v in d_value:
                        final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr_int),str(bins[x]), str(v[0]),str(v[1]),str(v[2]))
                        out.write(final)
                x += 1

write()
