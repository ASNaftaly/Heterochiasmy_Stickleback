#converting observed counts for 20 bins for ATAC data into proportions
#1. pull bin counts
#2. total bin counts
#3. calculate the total number of peaks observed across all chromosomes
#4. collapse bins to 10 bins (collapsed chromosome)
#5. calculate proportion of peaks in a given bin (0-0.5)
#4. write proportions to with only 10 bins (collapsed chromosome)
#to run script: python3 CalcProportion_PeakBins.py <observed bin counts in 20 bins> <output file with 10 bins in proportions not counts>
#Author: Alice Naftaly, August 2019

import sys


#pulls bins into a dictionary
#dictionary = key = bin number [0.05-1] and value = peak count for a given bin
#the input file has all autosomes (20 in total) so there should be 20 values for each dictionary key
def pull_peak_bins():
    bin_file = sys.argv[1]
    bin_dict = {}
    with open(bin_file, 'r') as bins:
        for line in bins:
            if line.startswith("Chr"):
                continue
            else:
                new_line = line.split()
                bin_number = new_line[1]
                peak_counts = new_line[4]
                if bin_number in bin_dict:
                    bin_dict[bin_number].append(peak_counts)
                elif bin_number not in bin_dict:
                    bin_dict.update({bin_number:[peak_counts]})
    return bin_dict

#calculates total peak counts across autosomes for a single bin
#returns dictionary with key = bin number and value = total count for bin
#these numbers will be the
def calc_total_peak_counts():
    bin_dict = pull_peak_bins()
    total_dict = {}
    for key in bin_dict:
        bin_number = key
        all_peak_counts = bin_dict[key]
        int_peak_counts = [int(i) for i in all_peak_counts]
        sum_counts = sum(int_peak_counts)
        total_dict.update({bin_number:sum_counts})
    return total_dict

#calculates the total number of peaks observed across all bins
#this is the denominator for the proportion
def total_observed_peaks():
    totals_dict = calc_total_peak_counts()
    total_observed_peaks = 0
    for key in totals_dict:
        single_bin_total = totals_dict[key]
        total_observed_peaks += single_bin_total
    return total_observed_peaks

#collapse bins from 20 bins to 10 bins
#returns list of peak count totals for each collapsed bins
def collapse_bins():
    bin_dict = calc_total_peak_counts()
    bin_1 = bin_dict["0.05"] + bin_dict["1.0"]
    bin_2 = bin_dict["0.1"] + bin_dict["0.95"]
    bin_3 = bin_dict["0.15"] + bin_dict["0.9"]
    bin_4 = bin_dict["0.2"] + bin_dict["0.85"]
    bin_5 = bin_dict["0.25"] + bin_dict["0.8"]
    bin_6 = bin_dict["0.3"] + bin_dict["0.75"]
    bin_7 = bin_dict["0.35"] + bin_dict["0.7"]
    bin_8 = bin_dict["0.4"] + bin_dict["0.65"]
    bin_9 = bin_dict["0.45"] + bin_dict["0.6"]
    bin_10 = bin_dict["0.5"] + bin_dict["0.55"]
    collapsed_bins = [bin_1, bin_2, bin_3, bin_4, bin_5, bin_6, bin_7, bin_8, bin_9, bin_10]
    return collapsed_bins

#calculates proportion of peaks in each bin (10 bins total)
#returns proportion as a percent
#returns list of proportions for 10 bins
def calc_proportions():
    total_count = total_observed_peaks()
    bin_counts = collapse_bins()
    bin_proportions = []
    for value in bin_counts:
        single_bin_proportion = round((value/total_count)*100,2)
        bin_proportions.append(single_bin_proportion)
    return bin_proportions

#writes output to file
#file format:
#Bin.Num    Peak.proportion
#only 10 bins; will not need to collapse bins in R or for making figure
def write():
    final_proportions = calc_proportions()
    bins = ("0.05", "0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5")
    output = sys.argv[2]
    with open(output, 'a') as out:
        header = "Bin.Num\tPeak.Proportion\n"
        out.write(header)
        for index, value in enumerate(bins):
            final = "%s\t%s\n" % (str(value),str(final_proportions[index]))
            out.write(final)

write()
