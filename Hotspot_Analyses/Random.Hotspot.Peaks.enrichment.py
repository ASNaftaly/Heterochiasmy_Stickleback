#comparing hotspots and peaks
#this is the random enrichment script that will be run 10,000 times
#to run script: python3 Random.Hotspot.Peaks.enrichment.py <peak positions from all gonads> <random hotspot file with chr \t spot start \t spot end> <permutation number> <output file: perm number \t num hotspots at TSSs>
#Author: Alice Naftaly, Oct 2020

import sys

#read in transcript positions
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

#read in random hotspots
#returns dictionary with key == chr number and value == [spot start, spot end]
def read_random_spots():
    spots_file = sys.argv[2]
    spots_dict = {}
    with open(spots_file, 'r') as spots:
        for line in spots:
            new_line = line.split()
            chr_num = new_line[0]
            spot_start = int(new_line[1])
            spot_end = int(new_line[2])
            dict_value = [spot_start, spot_end]
            if chr_num in spots_dict:
                spots_dict[chr_num].append(dict_value)
            elif chr_num not in spots_dict:
                spots_dict.update({chr_num:[dict_value]})
    return spots_dict

#count how many random hotspots overlap peaks
#returns total number of hotspots that overlap peaks
def compare_peaks_spots():
    peaks = read_peaks()
    hotspots = read_random_spots()
    hotspot_overlapping_peaks = []
    for chr in hotspots:
        single_chr_peaks = peaks[chr]
        single_chr_hotspots = hotspots[chr]
        for peak in single_chr_peaks
        :
            peak_start = peak[0]
            peak_end = peak[1]
            for ind, hot in enumerate(single_chr_hotspots):
                hot_start = hot[0]
                hot_end = hot[1]
                if hot_start <= peak_start and hot_start <= peak_end and hot_end >= peak_start and hot_end >= peak_end:
                    hotspot_overlapping_peaks.append(ind)
                elif hot_start > peak_start and hot_start < peak_end and hot_end > peak_start and hot_end > peak_end:
                    hotspot_overlapping_peaks.append(ind)
                elif hot_start < peak_start and hot_start < peak_end and hot_end > peak_start and hot_end < peak_end:
                    hotspot_overlapping_peaks.append(ind)
    set_hot_overlap = len(list(set(hotspot_overlapping_peaks)))
    return set_hot_overlap

#write output
def write():
    overlap = compare_peaks_spots()
    perm_num = sys.argv[3]
    output = sys.argv[4]
    with open(output, 'a') as out:
        final = "%s\t%s\n" % (str(perm_num), str(overlap))
        out.write(final)

write()
