#script to combine overlapping ATAC peaks between samples from the same tissues
#overlapping peaks will be defined as peaks that overlap by 1 or more bps
#will read in filtered xls file; format:
#chr  start   end     length  abs_summit      pileup  -log10(pvalue)  fold_enrichment -log10(qvalue)  name
#will need to compare chromosome by chromosome
#peaks can overlap in 4 ways if searching peak A against peaks from file B
#1. peak B overlaps the start of peak A
#2. peak B can be found inside of peak A
#3. peak B overlaps the end of peak A
#4. peak A can be found inside of peak B
#output formats:
#Unique peaks files: Chr \t Peak Start \t Peak End \n
#Summary Stats file: Chr \t New Peak Start (identifier) \t peak A size \t peak B size \t Final Peak Size\n
#combined peaks file: chr \t start \t end \t length \n
#to run script: python3 CombinedPeaks.py <peak file A (xls)> <peak file B (xls)> <Output: Unique Peaks from file A> <Output: Unique Peaks from file B> <Output: Summary stats file> <Output: Combined Peaks File>
#Author: Alice Naftaly, January 2022

import sys

#read in peak file A
#save positions in a dictionary with key == chr and value == [(peak start, peak end),...]
def read_file_A():
    A_peaks = {}
    a_file = sys.argv[1]
    with open(a_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("chr"):
                new_line = line.split("\t")
                if new_line[1] == "start":
                    continue
                else:
                    chr = new_line[0]
                    peak_start = int(new_line[1])
                    peak_end = int(new_line[2])
                    dict_value = [peak_start, peak_end]
                    if chr in A_peaks:
                        A_peaks[chr].append(dict_value)
                    elif chr not in A_peaks:
                        A_peaks.update({chr:[dict_value]})
    return A_peaks


#read in peak file B
#save positions in a dictionary with key == chr and value == [(peak start, peak end),...]
def read_file_B():
    B_peaks = {}
    b_file = sys.argv[2]
    with open(b_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("chr"):
                new_line = line.split("\t")
                if new_line[1] == "start":
                    continue
                else:
                    chr = new_line[0]
                    peak_start = int(new_line[1])
                    peak_end = int(new_line[2])
                    dict_value = [peak_start, peak_end]
                    if chr in B_peaks:
                        B_peaks[chr].append(dict_value)
                    elif chr not in B_peaks:
                        B_peaks.update({chr:[dict_value]})
    return B_peaks


#combine peaks that overlap based on the criteria above
#returns dictionary with overlapping peaks where key == chr and value == [[peak a start, peak a end, peak b start, peak b end],...]
#note one peak from file A can cover multiple peaks in file B depending on positioning and size of the peak. This will show up as multiple instances of a peak in the overlapping_peaks dictionary. I will need to condense these in another function. This may not be an issue with the filtered peaks.
def overlapping_peaks():
    a_peaks = read_file_A()
    b_peaks = read_file_B()
    overlapping_peaks = {}
    a_peaks_unique = {}
    b_peaks_unique = {}
    for chr in a_peaks:
        single_chr_a = a_peaks[chr]
        single_chr_b = b_peaks[chr]
        a_peaks_overlapping = []
        b_peaks_overlapping = []
        for peak_a in single_chr_a:
            peak_a_start = peak_a[0]
            peak_a_end = peak_a[1]
            for peak_b in single_chr_b:
                peak_b_start = peak_b[0]
                peak_b_end = peak_b[1]
                #overlap scenario 1:
                if peak_a_start > peak_b_start and peak_a_start <= peak_b_end and peak_a_end > peak_b_start and peak_a_end > peak_b_end:
                    a_peaks_overlapping.append(peak_a)
                    b_peaks_overlapping.append(peak_b)
                    dict_value = [peak_a_start, peak_a_end, peak_b_start, peak_b_end]
                    if chr in overlapping_peaks:
                        overlapping_peaks[chr].append(dict_value)
                    elif chr not in overlapping_peaks:
                        overlapping_peaks.update({chr:[dict_value]})
                #overlap scenario 2:
                elif peak_a_start <= peak_b_start and peak_a_start < peak_b_end and peak_a_end > peak_b_start and peak_a_end >= peak_b_end:
                    a_peaks_overlapping.append(peak_a)
                    b_peaks_overlapping.append(peak_b)
                    dict_value = [peak_a_start, peak_a_end, peak_b_start, peak_b_end]
                    if chr in overlapping_peaks:
                        overlapping_peaks[chr].append(dict_value)
                    elif chr not in overlapping_peaks:
                        overlapping_peaks.update({chr:[dict_value]})
                #overlap scenario 3:
                elif peak_a_start < peak_b_start and peak_a_start < peak_b_end and peak_a_end >= peak_b_start and peak_a_end < peak_b_end:
                    a_peaks_overlapping.append(peak_a)
                    b_peaks_overlapping.append(peak_b)
                    dict_value = [peak_a_start, peak_a_end, peak_b_start, peak_b_end]
                    if chr in overlapping_peaks:
                        overlapping_peaks[chr].append(dict_value)
                    elif chr not in overlapping_peaks:
                        overlapping_peaks.update({chr:[dict_value]})
                #overlap scenario 4:
                elif peak_a_start >= peak_b_start and peak_a_start < peak_b_end and peak_a_end > peak_b_start and peak_a_end <= peak_b_end:
                    a_peaks_overlapping.append(peak_a)
                    b_peaks_overlapping.append(peak_b)
                    dict_value = [peak_a_start, peak_a_end, peak_b_start, peak_b_end]
                    if chr in overlapping_peaks:
                        overlapping_peaks[chr].append(dict_value)
                    elif chr not in overlapping_peaks:
                        overlapping_peaks.update({chr:[dict_value]})
        for peak in single_chr_a:
            if peak in a_peaks_overlapping:
                continue
            else:
                if chr in a_peaks_unique:
                    a_peaks_unique[chr].append(peak)
                elif chr not in a_peaks_unique:
                    a_peaks_unique.update({chr:[peak]})
        for peak2 in single_chr_b:
            if peak2 in b_peaks_overlapping:
                continue
            else:
                if chr in b_peaks_unique:
                    b_peaks_unique[chr].append(peak2)
                elif chr not in b_peaks_unique:
                    b_peaks_unique.update({chr:[peak2]})
    #code to confirm all peaks are being sorted properly
    print(len(overlapping_peaks["chrVI"]))
    print(len(a_peaks["chrVI"]))
    print(len(b_peaks["chrVI"]))
    print(len(a_peaks_unique["chrVI"]))
    print(len(b_peaks_unique["chrVI"]))
    overlapping_peak_counts = 0
    unique_a_peaks = 0
    peak_A_start = []
    peak_B_start = []
    for k in overlapping_peaks["chrVI"]:
        peak_A_start.append(str(k[0]))
        peak_B_start.append(str(k[2]))
    #change the following lines to > 1 or > 2 etc to get full counts
    dups_A = {tuple(x) for x in peak_A_start if peak_A_start.count(x) > 2}
    print(len(dups_A))
    dups_B = {tuple(x) for x in peak_B_start if peak_B_start.count(x) > 2}
    print(len(dups_B))
    '''return overlapping_peaks, a_peaks_unique, b_peaks_unique'''

overlapping_peaks()
#write out unique peaks to separate file to go through as needed
def write_A_unique_peaks():
    overlapping, a_unique_peaks, b_unique_peaks = overlapping_peaks()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Chr\tPeak.Start\tPeak.End\n"
        out.write(header)
        for chr in a_unique_peaks:
            single_chr = a_unique_peaks[chr]
            for peak in single_chr:
                final = "%s\t%s\t%s\n" % (chr, str(peak[0]), str(peak[1]))
                out.write(final)

def write_B_unique_peaks():
    overlapping, a_unique_peaks, b_unique_peaks = overlapping_peaks()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Chr\tPeak.Start\tPeak.End\n"
        out.write(header)
        for chr in b_unique_peaks:
            single_chr = b_unique_peaks[chr]
            for peak in single_chr:
                final = "%s\t%s\t%s\n" % (chr, str(peak[0]), str(peak[1]))
                out.write(final)


#For overlapping peaks, I need to condense them into a single peak for analysis. Overlapping peaks should already be together in a list, but I'll need to combine the peaks that overlap multiple peaks. Will do this with the finalized peaks
#will return 2 dictionaries:
#summary stats = sizes of original peaks with new peak values
#combined peaks = has information for final xls form to build
def combine_peaks():
    overlapping, a_unique_peaks, b_unique_peaks = overlapping_peaks()
    summary_stats = {}
    combined_peaks = {}
    for chr in overlapping:
        single_chr = overlapping[chr]
        for p in single_chr:
            a_peak_start = p[0]
            a_peak_end = p[1]
            b_peak_start = p[2]
            b_peak_end = p[3]
            min_peak_start = min(p)
            max_peak_end = max(p)
            original_peak_sizes = [str(a_peak_end-a_peak_start), str(b_peak_end-b_peak_start)]
            final_peak_size = str(max_peak_end-min_peak_start)
            summary_value = [min_peak_start, original_peak_sizes, final_peak_size]
            combined_value = [str(min_peak_start), str(max_peak_end), final_peak_size]
            if chr in summary_stats:
                summary_stats[chr].append(summary_value)
            elif chr not in summary_stats:
                summary_stats.update({chr:[summary_value]})
            if chr in combined_peaks:
                combined_peaks[chr].append(combined_value)
            elif chr not in combined_peaks:
                combined_peaks.update({chr:[combined_value]})
    return summary_stats, combined_peaks

#write output for summary stats for records and final xls (tab delimited)

#summary stats:
def write_summary():
    summary_stats, combined_peaks = combine_peaks()
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = "Chr\tNew.Peak\tSample1.Peak.Size\tSample2.Peak.Size\tFinal.Peak.Size\n"
        out.write(header)
        for chr in summary_stats:
            single_chr = summary_stats[chr]
            for a in single_chr:
                new_peak = a[0]
                peak_sizes = a[1]
                final_peak_size = a[2]
                final = "%s\t%s\t%s\t%s\n" % (chr, new_peak, peak_sizes[0], peak_sizes[1],final_peak_size)
                out.write(final)

#combined peaks
def write_combined_peaks():
    summary_stats, combined_peaks = combine_peaks()
    output = sys.argv[6]
    with open(output, 'a') as out:
        header = "chr\tstart\tend\tlength\n"
        out.write(header)
        for chr in combined_peaks:
            single_chr = combined_peaks[chr]
            for b in single_chr:
                final = "%s\t%s\t%s\t%s\n" % (chr, b[0], b[1], b[2])
                out.write(final)


#run all functions:
def run_all():
    sample_1_unique_peaks = write_A_unique_peaks()
    sample_2_unique_peaks = write_B_unique_peaks()
    summary_stats = write_summary()
    combined_peaks = write_combined_peaks()

#run_all()
