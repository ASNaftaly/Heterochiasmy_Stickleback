#creating 20 bins across chromosome with 5 bins of each size on each size of the centromere for ATAC-seq (read coverage and peaks)
#this script handles peaks
#will first pull chromosome and centromere position and figure out size of bins for each arm
#then will pull number of peaks in each bin
#to run script: python3 Make20Bins_ATAC_peaks.py <centromere positions file with roman numerals> <chr number as roman numeral> <peaks file> <readcov file> <output file>
#Author: Alice Shanfelter, May 2019, edited September 2020

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
    peaks_list = []
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            chr_num_full = chr_num_roman + "\t"
            if line.startswith(chr_num_full):
                new_line = line.split()
                peak_start = new_line[1]
                peak_end = new_line[2]
                peaks_list.append(peak_start)
                peaks_list.append(peak_end)
    return peaks_list

#pulls length of chromosome from readcov file
#returns integer of chromosome length
def pull_chr_length():
    chr_length_file = sys.argv[4]
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

def count_peaks_arm1():
    bins = create_bins()
    arm1_bin_size = bins[0]
    peaks = pull_peaks()
    peaks_start = peaks[0::2]
    peaks_end = peaks[1::2]
    bin_dict = {}
    x = 0
    bin_start = 0
    bin_end = arm1_bin_size
    while x < 10:
        bin_count = 0
        for index, value in enumerate(peaks_start):
            if int(value) >= bin_start and int(value) <= bin_end and int(peaks_end[index]) <= bin_end:
                bin_count += 1
                #print("Peak in Bin")
                #print(value)
                #print(peaks_end[index])
            elif int(value) >= bin_start and int(value) <= bin_end and int(peaks_end[index]) >= bin_end:
                bin_count += 1
                #print("Peak overlap Bin")
                #print(value)
                #print(peaks_end[index])
            else:
                continue
        final = [bin_start,bin_end, bin_count]
        bin_dict.update({x:final})
        bin_start += arm1_bin_size
        bin_end += arm1_bin_size
        bin_count = 0
        x += 1
    return bin_dict

def count_peaks_arm2():
    bins = create_bins()
    arm2_bin_size = bins[1]
    peaks = pull_peaks()
    peaks_start = peaks[0::2]
    peaks_end = peaks[1::2]
    centromere_pos = pull_centromere()
    bin_dict = {}
    x = 10
    bin_start = centromere_pos
    bin_end = arm2_bin_size + centromere_pos
    while x < 20:
        bin_count = 0
        for index, value in enumerate(peaks_start):
            if int(value) >= bin_start and int(value) <= bin_end and int(peaks_end[index]) <= bin_end:
                bin_count += 1
                #print("Peak in Bin")
                #print(value)
                #print(peaks_end[index])
            elif int(value) >= bin_start and int(value) <= bin_end and int(peaks_end[index]) >= bin_end:
                bin_count += 1
                #print("Peak overlap Bin")
                #print(value)
                #print(peaks_end[index])
            else:
                continue
        final = [bin_start,bin_end, bin_count]
        bin_dict.update({x:final})
        bin_start += arm2_bin_size
        bin_end += arm2_bin_size
        bin_count = 0
        x += 1
    return bin_dict

#writing to one output file
def write():
    arm1_dict = count_peaks_arm1()
    arm2_dict = count_peaks_arm2()
    chr_num = sys.argv[2]
    bins = ["0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1.0"]
    output = sys.argv[5]
    with open(output,'a') as out:
        header = "Chr.Num\tBin.Num\tBin.Start\tBin.End\tNumber.Peaks\n"
        out.write(header)
        x = 0
        while x < 20:
            if x in arm1_dict:
                single_dict_value = arm1_dict[x]
                final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr_num),str(bins[x]), str(single_dict_value[0]),str(single_dict_value[1]),str(single_dict_value[2]))
                out.write(final)
                x += 1
            elif x in arm2_dict:
                single_dict_value = arm2_dict[x]
                final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr_num),str(bins[x]), str(single_dict_value[0]),str(single_dict_value[1]),str(single_dict_value[2]))
                out.write(final)
                x += 1

write()
