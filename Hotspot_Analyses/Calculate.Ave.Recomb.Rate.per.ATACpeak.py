#script to calculate average recombination rate across ATAC peaks
#recombination rates from hotspot paper and converted positions to V5 assembly
#Recombination rates file format: (bed file) chr \t start position \t end position \t recombination rate of region
#ATAC file format: (xls file) chr \t start pos \t end pos + other columsn with other data
#just need to compare positions and pull average recombination rate
#final file format: chr \t peak.start.position \t peak.end.position \t ave.recombination.rate.of.peak
#to run script: python3 Calculate.Ave.Recomb.Rate.per.ATACpeak.py <recomb rates bed file> <ATAC xls file> <output file>
#Author: Alice Naftaly, May 2021

import sys

#read in recombination rates file
#returns dictionary with key == chromosome number and value == list of [start position of region, end position of region, average recombination rate for region]
def read_rates():
    rates_file = sys.argv[1]
    rates_dict = {}
    with open(rates_file, 'r') as rates:
        for line in rates:
            new_line = line.split()
            chr_num = new_line[0]
            start_pos = int(new_line[1])
            end_pos = int(new_line[2])
            recomb_rate = float(new_line[3])
            dict_value = [start_pos, end_pos, recomb_rate]
            if chr_num in rates_dict:
                rates_dict[chr_num].append(dict_value)
            elif chr_num not in rates_dict:
                rates_dict.update({chr_num:[dict_value]})
    return rates_dict

#read in ATAC peaks
#returns dictionary with key == chromosome number and value == list of peaks [start pos, end pos]
def read_ATAC():
    atac_file = sys.argv[2]
    atac_dict = {}
    with open(atac_file, 'r') as atac:
        for line in atac:
            new_line = line.split()
            chr_num = new_line[0]
            peak_start = int(new_line[1])
            peak_end = int(new_line[2])
            dict_value = [peak_start, peak_end]
            if chr_num in atac_dict:
                atac_dict[chr_num].append(dict_value)
            elif chr_num not in atac_dict:
                atac_dict.update({chr_num:[dict_value]})
    return atac_dict

##calculating average recombination rate per peak
def calc_ave_rate():
    atac_peaks = read_ATAC()
    recomb_rates = read_rates()
    peak_rates_dict = {}
    for chr in atac_peaks:
        single_chr_peaks = atac_peaks[chr]
        single_chr_rates = recomb_rates[chr]
        single_chr_combined = {}
        for peak in single_chr_peaks:
            peak_start = peak[0]
            peak_end = peak[1]
            for region in single_chr_rates:
                region_start = region[0]
                region_end = region[1]
                region_recomb_rate = region[2]
                #records peaks that fall fully within a single region
                if peak_start >= region_start and peak_start <= region_end and peak_end >= region_start and peak_end <= region_end:
                    value = [peak_start, peak_end, region_recomb_rate]
                    if peak_start in single_chr_combined:
                        single_chr_combined[peak_start].append(value)
                    elif peak_start not in single_chr_combined:
                        single_chr_combined.update({peak_start:[value]})
                #records when peak partially overlaps region (front half)
                elif peak_start >= region_start and peak_start <= region_end and peak_end >= region_start and peak_end >= region_end:
                    value = [peak_start, peak_end, region_recomb_rate, region_end - peak_start]
                    if peak_start in single_chr_combined:
                        single_chr_combined[peak_start].append(value)
                    elif peak_start not in single_chr_combined:
                        single_chr_combined.update({peak_start:[value]})
                #records when peak partially overlaps region (end half)
                elif peak_start <= region_start and peak_start <= region_end and peak_end >= region_start and peak_end <= region_end:
                    value = [peak_start, peak_end, region_recomb_rate, peak_end - region_start]
                    if peak_start in single_chr_combined:
                        single_chr_combined[peak_start].append(value)
                    elif peak_start not in single_chr_combined:
                        single_chr_combined.update({peak_start:[value]})
        for key in single_chr_combined:
            single_peak = single_chr_combined[key]
            if len(single_peak) == 1:
                if chr in peak_rates_dict:
                    peak_rates_dict[chr].append(single_peak[0])
                elif chr not in peak_rates_dict:
                    peak_rates_dict.update({chr:[single_peak[0]]})
            elif len(single_peak) == 2:
                peak_region_1_overlap = single_peak[0]
                peak_region_2_overlap = single_peak[1]
                average_rate = ((peak_region_1_overlap[2]*peak_region_1_overlap[3]) + (peak_region_2_overlap[2]*peak_region_2_overlap[3]))/(peak_region_1_overlap[1]-peak_region_1_overlap[0])
                final_value = [peak_region_1_overlap[0], peak_region_1_overlap[1], average_rate]
                if chr in peak_rates_dict:
                    peak_rates_dict[chr].append(final_value)
                elif chr not in peak_rates_dict:
                    peak_rates_dict.update({chr:[final_value]})
    return peak_rates_dict

#write output
def write():
    output = sys.argv[3]
    combined_dict = calc_ave_rate()
    with open(output, 'a') as out:
        header = "Chr\tPeak.Start\tPeak.End\tAve.Recomb.Rate\n"
        out.write(header)
        for key in combined_dict:
            single_chr = combined_dict[key]
            for peak in single_chr:
                final = "%s\t%s\t%s\t%s\n" % (key, str(peak[0]), str(peak[1]), str(peak[2]))
                out.write(final)


write()
