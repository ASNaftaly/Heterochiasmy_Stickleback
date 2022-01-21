#removing duplicate peaks from filtered combined peaks files
#will first condense duplicate peaks from file A and then duplicate peaks from file B
#where duplicate peaks while be condensed into a single peak
#to run script: python3 CondenseDuplicatePeaks.py <filtered combined peaks file> <Output: final filtered combined peaks file>
#Author: Alice Naftaly, Jan 2022

import sys

#read in filtered peaks file
#format: chr \t start \t end \t length
def read_peaks():
    peaks_file = sys.argv[1]
    peak_dict = {}
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            new_line = line.split("\t")
            if new_line[1] == "start":
                continue
            else:
                chr = new_line[0]
                peak_start = int(new_line[1])
                peak_end = int(new_line[2])
                dict_value = [peak_start, peak_end]
                if chr in peak_dict:
                    peak_dict[chr].append(dict_value)
                elif chr not in peak_dict:
                    peak_dict.update({chr:[dict_value]})
    return peak_dict


#condense peak duplicates based on start position
#returns a dictionary with key == chr and value == [[peak start, peak end],...]
def condense_start_pos():
    peaks = read_peaks()
    condensed_start_peaks = {}
    for chr in peaks:
        single_chr = peaks[chr]
        start_dict = {}
        for v in single_chr:
            start = v[0]
            end = v[1]
            start_value = [start, end]
            if start in start_dict:
                start_dict[start].append(start_value)
            elif start not in start_dict:
                start_dict.update({start:[start_value]})
        for key in start_dict:
            single_key = start_dict[key]
            if len(single_key) == 1:
                peak_start = single_key[0][0]
                peak_end = single_key[0][1]
                dict_value = [peak_start, peak_end]
                if chr in condensed_start_peaks:
                    condensed_start_peaks[chr].append(dict_value)
                elif chr not in condensed_start_peaks:
                    condensed_start_peaks.update({chr:[dict_value]})
            elif len(single_key) > 1:
                start_list = [item[0] for item in single_key]
                end_list = [item[1] for item in single_key]
                final_start = min(start_list)
                final_end = max(end_list)
                dict_value = [final_start, final_end]
                if chr in condensed_start_peaks:
                    condensed_start_peaks[chr].append(dict_value)
                elif chr not in condensed_start_peaks:
                    condensed_start_peaks.update({chr:[dict_value]})
    return condensed_start_peaks


#condense peak duplicates based on end position
#returns a dictionary with key == chr and value == [[peak start, peak end],...]
def condense_end_pos():
    peaks = condense_start_pos()
    condensed_end_peaks = {}
    for chr in peaks:
        single_chr = peaks[chr]
        end_dict = {}
        for v in single_chr:
            start = v[0]
            end = v[1]
            end_value = [start, end]
            if end in end_dict:
                end_dict[end].append(end_value)
            elif end not in end_dict:
                end_dict.update({end:[end_value]})
        for key in end_dict:
            single_key = end_dict[key]
            if len(single_key) == 1:
                peak_start = single_key[0][0]
                peak_end = single_key[0][1]
                dict_value = [peak_start, peak_end]
                if chr in condensed_end_peaks:
                    condensed_end_peaks[chr].append(dict_value)
                elif chr not in condensed_end_peaks:
                    condensed_end_peaks.update({chr:[dict_value]})
            elif len(single_key) > 1:
                start_list = [item[0] for item in single_key]
                end_list = [item[1] for item in single_key]
                final_start = min(start_list)
                final_end = max(end_list)
                dict_value = [final_start, final_end]
                if chr in condensed_end_peaks:
                    condensed_end_peaks[chr].append(dict_value)
                elif chr not in condensed_end_peaks:
                    condensed_end_peaks.update({chr:[dict_value]})
    return condensed_end_peaks


#write final condensed peak file
def write():
    peaks = condense_end_pos()
    output = sys.argv[2]
    with open(output, 'a') as out:
        header = "chr\tstart\tend\tlength\n"
        for chr in peaks:
            single_chr = peaks[chr]
            for peak in single_chr:
                p_start = peak[0]
                p_end = peak[1]
                p_length = p_end - p_start
                final = "%s\t%s\t%s\t%s\n" % (chr, str(p_start), str(p_end), str(p_length))
                out.write(final)


write()
