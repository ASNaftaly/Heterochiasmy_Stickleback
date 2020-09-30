#Calculation proportions for 10k random permutations to create range of values expected if peaks were randomly distributed
#permutation summary files format: each line has bin number and then 10,000 total count values across autosomes
#1. read permutation summary file and observed 20 bins counts file
#2. calculate total peaks observed from observed counts file
#3. collapse bins from 20 bins to 10 bins
#4. calculate proportion for each permutation for each bin
#5. write proportions to new file with only 10 bins
#to run script: python3 CalcProportion_random_permutations_atac.py <permutations summary file> <observed bin counts for 20 bins> <output file for 10 bins with 10k permutations with proportions>
#Author: Alice Naftaly, August 2019

import sys

def pull_permutations():
    permutations_file = sys.argv[1]
    permutations_dict = {}
    with open(permutations_file, 'r') as permutations:
        for line in permutations:
            new_line = line.split()
            bin_number = new_line[0]
            permutation_counts = new_line[1:len(new_line)]
            permutations_dict.update({bin_number:permutation_counts})
    return permutations_dict

#pulls bins into a dictionary
#dictionary = key = bin number [0.05-1] and value = peak count for a given bin
#the input file has all autosomes (20 in total) so there should be 20 values for each dictionary key
def pull_peak_bins():
    bin_file = sys.argv[2]
    bin_dict = {}
    with open(bin_file, 'r') as bins:
        for line in bins:
            if line.startswith("Bin.Num"):
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

#calculates the total number of peaks observed across all bins
#this is the denominator for the proportion
def total_observed_peaks():
    totals_dict = pull_peak_bins()
    total_observed_peaks = 0
    for key in totals_dict:
        single_bin = totals_dict[key]
        int_single_bin = [int(i) for i in single_bin]
        sum_single_bin = sum(int_single_bin)
        total_observed_peaks += sum_single_bin
    return total_observed_peaks


#collapse permutations from 20 bins to 10 bins
#this function is a bit longer and has more hard coding than I would like
#might adjust this later
#returns dictionary with 10 bins where key = bin number (1-10) and value = 10k permutations for total peak counts for collapsed bin
def collapse_bins():
    permutations = pull_permutations()
    bin_1 = permutations["0.05"]
    bin_2 = permutations["0.10"]
    bin_3 = permutations["0.15"]
    bin_4 = permutations["0.20"]
    bin_5 = permutations["0.25"]
    bin_6 = permutations["0.30"]
    bin_7 = permutations["0.35"]
    bin_8 = permutations["0.40"]
    bin_9 = permutations["0.45"]
    bin_10 = permutations["0.50"]
    bin_11 = permutations["0.55"]
    bin_12 = permutations["0.60"]
    bin_13 = permutations["0.65"]
    bin_14 = permutations["0.70"]
    bin_15 = permutations["0.75"]
    bin_16 = permutations["0.80"]
    bin_17 = permutations["0.85"]
    bin_18 = permutations["0.90"]
    bin_19 = permutations["0.95"]
    bin_20 = permutations["1.0"]
    collapsed_bin_1 = []
    collapsed_bin_2 = []
    collapsed_bin_3 = []
    collapsed_bin_4 = []
    collapsed_bin_5 = []
    collapsed_bin_6 = []
    collapsed_bin_7 = []
    collapsed_bin_8 = []
    collapsed_bin_9 = []
    collapsed_bin_10 = []
    collapsed_dict = {}
    x = 0
    while x < 10000:
        collapsed_bin_1.append(int(bin_1[x]) + int(bin_20[x]))
        collapsed_bin_2.append(int(bin_2[x]) + int(bin_19[x]))
        collapsed_bin_3.append(int(bin_3[x]) + int(bin_18[x]))
        collapsed_bin_4.append(int(bin_4[x]) + int(bin_17[x]))
        collapsed_bin_5.append(int(bin_5[x]) + int(bin_16[x]))
        collapsed_bin_6.append(int(bin_6[x]) + int(bin_15[x]))
        collapsed_bin_7.append(int(bin_7[x]) + int(bin_14[x]))
        collapsed_bin_8.append(int(bin_8[x]) + int(bin_13[x]))
        collapsed_bin_9.append(int(bin_9[x]) + int(bin_12[x]))
        collapsed_bin_10.append(int(bin_10[x]) + int(bin_11[x]))
        x +=1
    collapsed_dict.update({"1":collapsed_bin_1})
    collapsed_dict.update({"2":collapsed_bin_2})
    collapsed_dict.update({"3":collapsed_bin_3})
    collapsed_dict.update({"4":collapsed_bin_4})
    collapsed_dict.update({"5":collapsed_bin_5})
    collapsed_dict.update({"6":collapsed_bin_6})
    collapsed_dict.update({"7":collapsed_bin_7})
    collapsed_dict.update({"8":collapsed_bin_8})
    collapsed_dict.update({"9":collapsed_bin_9})
    collapsed_dict.update({"10":collapsed_bin_10})
    return collapsed_dict

#calculations proportions of peaks in each bin
#returns dictionary with key = bin number (1-10) and value = 10k permutations in proportions
def calc_proportions():
    permutations_dict = collapse_bins()
    total_count = total_observed_peaks()
    proportions_dict = {}
    for key in permutations_dict:
        single_key = permutations_dict[key]
        proportions = [round((int(i)/total_count)*100,2) for i in single_key]
        proportions_dict.update({key:proportions})
    return proportions_dict


#writes output file
#final output has proportions of peaks in each bin for 10k permutations
#output format: Bin Number 10k permutation proportions separated by tabs
def write():
    proportions = calc_proportions()
    bins = ("0.05", "0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5")
    output = sys.argv[3]
    with open(output, 'a') as out:
        for index, value in enumerate(bins):
            single_bin = proportions[str(index+1)]
            str_single_bin = [str(i) for i in single_bin]
            join_values = "\t".join(str_single_bin)
            final = str(value) + "\t" + join_values + "\n"
            out.write(final)

write()
