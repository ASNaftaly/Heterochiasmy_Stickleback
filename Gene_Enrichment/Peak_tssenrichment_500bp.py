#Checking peaks identified in Liver and Testis
#First pull peaks associated within 500bp of TSSs, then write peak position
#will first run analysis without any filtering outside pipeline
#returns list of ensembl transcript ids that are enriched around peaks for a specific chromosome
#to run script: python3 Peak_tssenrichment_500bp.py <tss file> <peaks file> <chr num as roman numeral> <output file>
#Author: Alice Shanfelter, March 2019, edited August 2019

import sys
import time

#pulling TSSs from psl format (BLAT ENSEMBl97 to GLAZER)
#returns list of transcription start sites
def pull_TSS():
    tss_file = sys.argv[1]
    chr_num = sys.argv[3]
    transcript_starts = []
    with open(tss_file,'r') as tss:
        for line in tss:
            new_line = line.split()
            chrnum = new_line[13]
            if chrnum == chr_num:
                transcript_starts.append(new_line[15])
    return transcript_starts

#pull peak start and end positions
#returns dictionary with key = peak start position and value = peak end position
def pull_peaks():
    peak_file = sys.argv[2]
    chr_num = sys.argv[3]
    chrnum= chr_num + "\t"
    peak_dict = {}
    with open(peak_file,'r') as peaks:
        for line in peaks:
            if line.startswith(chrnum):
                new_line = line.split()
                start_pos = new_line[1]
                end_pos = new_line[2]
                peak_dict.update({start_pos:end_pos})
    return peak_dict

#compares peaks to TSSs
#if peaks is within 500bp of TSS, this will be counted as overlapping
def compare_peaks():
    tss = pull_TSS()
    peaks = pull_peaks()
    peak_start = []
    peak_end = []
    enriched_peaks = []
    for key in peaks:
        peak_start.append(key)
        peak_end.append(peaks[key])
    x = 0
    while x < len(peak_start):
        for v in tss:
            upstream_tss = int(v) - 500
            downstream_tss = int(v) + 500
            #if peak overlaps upstream region of TSS, but doesn't extend past the downstream TSS region
            if upstream_tss <= int(peak_start[x]) and upstream_tss <= int(peak_end[x]) and downstream_tss >= int(peak_start[x]) and downstream_tss <= int(peak_end[x]):
                enriched_peaks.append(peak_start[x])
                enriched_peaks.append(peak_end[x])
                enriched_peaks.append(v)
            elif upstream_tss >= int(peak_start[x]) and upstream_tss <= int(peak_end[x]) and downstream_tss >= int(peak_start[x]) and downstream_tss <= int(peak_end[x]):
                enriched_peaks.append(peak_start[x])
                enriched_peaks.append(peak_end[x])
                enriched_peaks.append(v)
            elif upstream_tss >= int(peak_start[x]) and upstream_tss <= int(peak_end[x]) and downstream_tss >= int(peak_start[x]) and downstream_tss >= int(peak_end[x]):
                enriched_peaks.append(peak_start[x])
                enriched_peaks.append(peak_end[x])
                enriched_peaks.append(v)
        x += 1
    return enriched_peaks

#writes peak positions for those located around TSSs
def write_peaks():
    peaks = compare_peaks()
    peak_start = peaks[0::3]
    peak_end = peaks[1::3]
    tss_pos = peaks[2::3]
    header = "Peak.Start" + "\t" + "Peak.End" + "\t" + "TSS.pos" + "\n"
    output = sys.argv[4]
    with open(output, 'a') as out:
        out.write(header)
        for index, value in enumerate(peak_start):
            final = "%s\t%s\t%s\n" % (value,peak_end[index],tss_pos[index])
            out.write(final)
write_peaks()
