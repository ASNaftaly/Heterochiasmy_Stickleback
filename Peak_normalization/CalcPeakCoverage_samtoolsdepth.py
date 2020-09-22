#averaging bp coverage after samtools depth function is used for peaks
#needd to read in peak bed file and samtools depth output
#will produced bedgraph format as output
#chr  peak_start  peak_end  average_coverage
#to run script: python3 CalcPeakCoverage_samtoolsdepth.py <peak bed file> <samtools depth output> <chr_num> <output file>
#Author: Alice Naftaly, Sept 2020

import sys

#read in peaks bed file
#returns dictionary with key == chr num and value == [peak start, peak end, peak identifier]
def read_peaks():
    peaks_file = sys.argv[1]
    peaks_dict = {}
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            new_line = line.split()
            chr_num = new_line[0]
            peak_start = int(new_line[1])
            peak_end = int(new_line[2])
            peak_iden = new_line[3]
            dict_value = [peak_start, peak_end, peak_iden]
            if chr_num in peaks_dict:
                peaks_dict[chr_num].append(dict_value)
            elif chr_num not in peaks_dict:
                peaks_dict.update({chr_num:[dict_value]})
    return peaks_dict

#read in samtools output
def read_depth():
    depth_file = sys.argv[2]
    depth_dict = {}
    with open(depth_file, 'r') as depth:
        for line in depth:
            new_line = line.split()
            chr_num = new_line[0]
            position = int(new_line[1])
            coverage = int(new_line[2])
            dict_value = [position, coverage]
            if chr_num in depth_dict:
                depth_dict[chr_num].append(dict_value)
            elif chr_num not in depth_dict:
                depth_dict.update({chr_num:[dict_value]})
    return depth_dict


#average coverage across peaks:
def average_coverage():
    peaks = read_peaks()
    depth = read_depth()
    final_cov_list = []
    chr_num = sys.argv[3]
    single_chr_peaks = peaks[chr_num]
    single_chr_depth = depth[chr_num]
    for single_peak in single_chr_peaks:
        peak_start = single_peak[0]
        peak_end = single_peak[1]
        peak_iden = single_peak[2]
        peak_total_coverage = []
        for single_bp in single_chr_depth:
            bp_position = single_bp[0]
            coverage = single_bp[1]
            if peak_start <= bp_position <= peak_end:
                peak_total_coverage.append(coverage)
        ave_coverage = int(sum(peak_total_coverage)/len(peak_total_coverage))
        dict_value = [str(peak_start), str(peak_end), str(ave_coverage)]
        final_cov_list.append(dict_value)
    return final_cov_list


#write output as bedgraph format
def write():
    final_cov = average_coverage()
    chr_num = sys.argv[3]
    output = sys.argv[4]
    with open(output, 'a') as out:
        for value in final_cov:
            final = "%s\t%s\t%s\t%s\n" % (str(chr_num), value[0], value[1], value[2])
            out.write(final)


write()
