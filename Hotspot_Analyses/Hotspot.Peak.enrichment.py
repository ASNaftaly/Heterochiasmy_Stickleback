#comparing hotspots and peaks
#asking how many peaks fall at hotspots
#to run script: python3 Hotspot.Peak.enrichment.py <peak positions in xls file; where both ovaries and testes have been combined> <hotspots file in bed format; converted to V5 positions> <output file; has one hotspot per line that has a peak> >> <standard output has total number of peaks overlapping hotspots>
#Author: Alice Naftaly, Oct 2020

import sys

#read in peak positions
#returns dictionary with key == chr number and value = [peak start, peak end] * all peaks for chromosome
def read_peaks():
    peaks_file = sys.argv[1]
    peaks_dict = {}
    with open(peaks_file ,'r') as peaks:
        for line in peaks:
            new_line = line.split()
            chr_num = new_line[0]
            peak_start = int(new_line[1])
            peak_end = int(new_line[2])
            dict_value = [peak_start, peak_end]
            if chr_num in peaks_dict:
                peaks_dict[chr_num].append(dict_value)
            elif chr_num not in peaks_dict:
                peaks_dict.update({chr_num:[dict_value]})
    return peaks_dict

#read in hotspots
#returns dictionary with key == chr number and value == [hotspot start, hotspot end, hotspots ID]
def read_hotspots():
    hotspots_file = sys.argv[2]
    hotspots_dict = {}
    with open(hotspots_file, 'r') as hotspots:
        for line in hotspots:
            new_line = line.split()
            chr_num = new_line[0]
            hot_start = int(new_line[1])
            hot_end = int(new_line[2])
            hot_id = new_line[3]
            dict_value = [hot_start, hot_end, hot_id]
            if chr_num in hotspots_dict:
                hotspots_dict[chr_num].append(dict_value)
            elif chr_num not in hotspots_dict:
                hotspots_dict.update({chr_num:[dict_value]})
    return hotspots_dict


#count how many hotspots have peaks
def compare_peaks_hotspots():
    peaks = read_peaks()
    hotspots = read_hotspots()
    hotspots_overlapping_peaks = []
    peak_overlap_total = []
    for chr in hotspots:
        single_chr_peaks = peaks[chr]
        single_chr_hotspots = hotspots[chr]
        for peak in single_chr_peaks:
            peak_start = peak[0]
            peak_end = peak[1]
            peak_overlap = 0
            for hot in single_chr_hotspots:
                hot_start = hot[0]
                hot_end = hot[1]
                hot_id = hot[2]
                if hot_start <= peak_start and hot_start <= peak_end and hot_end >= peak_start and hot_end >= peak_end:
                    peak_overlap += 1
                    hotspots_overlapping_peaks.append(hot_id)
                elif hot_start > peak_start and hot_start < peak_end and hot_end > peak_start and hot_end > peak_end:
                    peak_overlap += 1
                    hotspots_overlapping_peaks.append(hot_id)
                elif hot_start < peak_start and hot_start < peak_end and hot_end > peak_start and hot_end < peak_end:
                    peak_overlap += 1
                    hotspots_overlapping_peaks.append(hot_id)
            if peak_overlap > 0:
                peak_overlap_total.append(peak_overlap)
    sum_peak_overlap = sum(peak_overlap_total)
    set_hot_overlap = list(set(hotspots_overlapping_peaks))
    return sum_peak_overlap, set_hot_overlap

#write output
def write():
    peak_overlap, hot_overlap = compare_peaks_hotspots()
    final_peak_overlap = "Total Peaks that overlap hotspots: %s\n" % str(peak_overlap)
    print(final_peak_overlap)
    output = sys.argv[3]
    with open(output, 'a') as out:
        for hot in hot_overlap:
            final_hot = "%s\n" % str(hot)
            out.write(final_hot)

write()
