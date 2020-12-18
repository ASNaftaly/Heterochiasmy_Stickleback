#collapse 20 bins across chromosomes
#follows MakeBins_peaksnearTSS_atac.py
#to run script: python3 CollapseBins_PeaksnotatTSSs.py <peaks not at tss in 20 bins; output of MakeBins_peaksnotatTSSs.py> <output>
#Author: Alice Naftaly, Oct 2020

import sys

#read in tsss in bins
#returns dictionary with key = bin number (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0) and value = list of counts (number of tss in peaks for a given bin)
def read_peaks_in_bins():
    peaks_file = sys.argv[1]
    bin_dict = {}
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("Chr.Num"):
                continue
            else:
                new_line = line.split()
                bin_number = new_line[1]
                bin_counts = int(new_line[4])
                if bin_number in bin_dict:
                    bin_dict[bin_number].append(bin_counts)
                elif bin_number not in bin_dict:
                    bin_dict.update({bin_number:[bin_counts]})
    return bin_dict

#total counts per bin
#returns dictionary with key = bin number and value = total counts for bin
def total_bins():
    bin_dict = read_peaks_in_bins()
    bin_totals = {}
    for bin in bin_dict:
        bin_total = sum(bin_dict[bin])
        bin_totals.update({bin:bin_total})
    return bin_totals

#calculate total number of peaks in tss for calculating percentages
def calc_sample_total():
    bin_dict = read_peaks_in_bins()
    total_peaks = 0
    for bin in bin_dict:
        bin_total = sum(bin_dict[bin])
        total_peaks += bin_total
    return total_peaks

#collapse bins
#returns dictionary with key == bin number (1-10) and value == [bin total counts, bin percents]
def collapse_bins():
    bin_dict = total_bins()
    total_peaks = calc_sample_total()
    collapsed_bins_dict = {}
    combined_bin_1_counts = bin_dict["0.05"] + bin_dict["1.0"]
    bin_1_percent = round((combined_bin_1_counts/total_peaks)*100,2)
    bin_1 = [combined_bin_1_counts, bin_1_percent]
    combined_bin_2_counts = bin_dict["0.1"] + bin_dict["0.95"]
    bin_2_percent = round((combined_bin_2_counts/total_peaks)*100,2)
    bin_2 = [combined_bin_2_counts, bin_2_percent]
    combined_bin_3_counts = bin_dict["0.15"] + bin_dict["0.9"]
    bin_3_percent = round((combined_bin_3_counts/total_peaks)*100,2)
    bin_3 = [combined_bin_3_counts, bin_3_percent]
    combined_bin_4_counts = bin_dict["0.2"] + bin_dict["0.85"]
    bin_4_percent = round((combined_bin_4_counts/total_peaks)*100,2)
    bin_4 = [combined_bin_4_counts, bin_4_percent]
    combined_bin_5_counts = bin_dict["0.25"] + bin_dict["0.8"]
    bin_5_percent = round((combined_bin_5_counts/total_peaks)*100,2)
    bin_5 = [combined_bin_5_counts, bin_5_percent]
    combined_bin_6_counts = bin_dict["0.3"] + bin_dict["0.75"]
    bin_6_percent = round((combined_bin_6_counts/total_peaks)*100,2)
    bin_6 = [combined_bin_6_counts, bin_6_percent]
    combined_bin_7_counts = bin_dict["0.35"] + bin_dict["0.7"]
    bin_7_percent = round((combined_bin_7_counts/total_peaks)*100,2)
    bin_7 = [combined_bin_7_counts, bin_7_percent]
    combined_bin_8_counts = bin_dict["0.4"] + bin_dict["0.65"]
    bin_8_percent = round((combined_bin_8_counts/total_peaks)*100,2)
    bin_8 = [combined_bin_8_counts, bin_8_percent]
    combined_bin_9_counts = bin_dict["0.45"] + bin_dict["0.6"]
    bin_9_percent = round((combined_bin_9_counts/total_peaks)*100,2)
    bin_9 = [combined_bin_9_counts, bin_9_percent]
    combined_bin_10_counts = bin_dict["0.5"] + bin_dict["0.55"]
    bin_10_percent = round((combined_bin_10_counts/total_peaks)*100,2)
    bin_10 = [combined_bin_10_counts, bin_10_percent]
    collapsed_bins_dict.update({"1":bin_1})
    collapsed_bins_dict.update({"2":bin_2})
    collapsed_bins_dict.update({"3":bin_3})
    collapsed_bins_dict.update({"4":bin_4})
    collapsed_bins_dict.update({"5":bin_5})
    collapsed_bins_dict.update({"6":bin_6})
    collapsed_bins_dict.update({"7":bin_7})
    collapsed_bins_dict.update({"8":bin_8})
    collapsed_bins_dict.update({"9":bin_9})
    collapsed_bins_dict.update({"10":bin_10})
    return collapsed_bins_dict

#write output
#output format: Bin.Number (1-10) \t Total.Peaks.in.TSS.Counts \t Percent.Peaks.in.TSS \n
def write():
    bin_percents = collapse_bins()
    output = sys.argv[2]
    with open(output, 'a') as out:
        header = "Bin.Number\tTotal.Peak.Counts\tPercent\n"
        out.write(header)
        for bin in bin_percents:
            final = "%s\t%s\t%s\n" % (str(bin), str(bin_percents[bin][0]), str(bin_percents[bin][1]))
            out.write(final)

write()
