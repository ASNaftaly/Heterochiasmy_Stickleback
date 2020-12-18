#need to pull peaks that are not at TSSs (from Pull.Peaks.no.TSSs.py) and separate these into peaks at TEs and peaks not at TEs or TSSs
#overlap with TSSs is within 500bp of TSS
#TE overlap is from Identify.TEs.at.Peaks.py
#need to read in output from both of these scripts
#to run script: python3 Split.Peaks.at.TEs.vs.not.at.TEs.or.TSSs.py <output from Pull.Peaks.no.TSSs.py> <DNA TE overlap > <LINE TE overlap> <SINE TE overlap> <LTR TE overlap> <RC TE overlap> <output file>
#Author: Alice Naftaly, Sept 2020

import sys

#read in file from Pull.Peaks.no.TSSs.py
#format: chr num \t peak start \t peak end
#returns dictionary as key == chr number and value == [peak start, peak end]
def read_peaks_not_at_tss():
    tss_file = sys.argv[1]
    peaks_not_at_tss = {}
    with open(tss_file, 'r') as tss:
        for line in tss:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[1])
                peak_end = int(new_line[2])
                peak = [peak_start, peak_end]
                if chr_num in peaks_not_at_tss:
                    peaks_not_at_tss[chr_num].append(peak)
                elif chr_num not in peaks_not_at_tss:
                    peaks_not_at_tss.update({chr_num:[peak]})
    return peaks_not_at_tss



#read in TE peak overlap
#format: 2 headers (1 TE type overlapping sigpeaks and header with column names)
#chr num, Overlap.type, TE.name TE.start TE.end Peak.start peak.end
def read_te_overlap():
    dna_te_file = sys.argv[2]
    line_te_file = sys.argv[3]
    sine_te_file = sys.argv[4]
    ltr_te_file = sys.argv[5]
    rc_te_file = sys.argv[6]
    peaks_overlapping_with_tes = {}
    with open(dna_te_file, 'r') as dna_te, open(line_te_file, 'r') as line_te, open(sine_te_file, 'r') as sine_te, open(ltr_te_file, 'r') as ltr_te, open(rc_te_file, 'r') as rc_te:
        for line in dna_te:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[5])
                peak_end = int(new_line[6])
                peak = [peak_start, peak_end]
                if chr_num in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes[chr_num].append(peak)
                elif chr_num not in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes.update({chr_num:[peak]})
        for line in line_te:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[5])
                peak_end = int(new_line[6])
                peak = [peak_start, peak_end]
                if chr_num in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes[chr_num].append(peak)
                elif chr_num not in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes.update({chr_num:[peak]})
        for line in sine_te:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[5])
                peak_end = int(new_line[6])
                peak = [peak_start, peak_end]
                if chr_num in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes[chr_num].append(peak)
                elif chr_num not in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes.update({chr_num:[peak]})
        for line in ltr_te:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[5])
                peak_end = int(new_line[6])
                peak = [peak_start, peak_end]
                if chr_num in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes[chr_num].append(peak)
                elif chr_num not in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes.update({chr_num:[peak]})
        for line in rc_te:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[5])
                peak_end = int(new_line[6])
                peak = [peak_start, peak_end]
                if chr_num in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes[chr_num].append(peak)
                elif chr_num not in peaks_overlapping_with_tes:
                    peaks_overlapping_with_tes.update({chr_num:[peak]})
    return peaks_overlapping_with_tes


#remove peaks that overlap with TEs from peaks that do not overlap with TSSs to get peaks that do not overlap TEs or TSSs
#returns dictionary with key = chr number and value == list of peaks
def filter_peaks():
    peaks_overlapping_with_tes = read_te_overlap()
    peaks_not_overlapping_tss = read_peaks_not_at_tss()
    peaks_not_at_tsss_or_tes = {}
    for chr in peaks_not_overlapping_tss:
        single_peaks_not_at_tsss = peaks_not_overlapping_tss[chr]
        single_peaks_at_tes = peaks_overlapping_with_tes[chr]
        for peak in single_peaks_not_at_tsss:
            if peak in single_peaks_at_tes:
                continue
            elif peak not in single_peaks_at_tes:
                if chr in peaks_not_at_tsss_or_tes:
                    peaks_not_at_tsss_or_tes[chr].append(peak)
                elif chr not in peaks_not_at_tsss_or_tes:
                    peaks_not_at_tsss_or_tes.update({chr:[peak]})
    return peaks_not_at_tsss_or_tes

#write output
#format:
#chr num \t peak start \t peak end
def write():
    peaks_not_at_tss_or_te = filter_peaks()
    output = sys.argv[7]
    chrs = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrXVII", 'chrXVIII', "chrXX", "chrXXI"]
    with open(output, 'a') as out:
        for chr in chrs:
            peak_list = []
            single_chr = peaks_not_at_tss_or_te[chr]
            for peak in single_chr:
                if peak in peak_list:
                    continue
                elif peak not in peak_list:
                    final = "%s\t%s\t%s\n" % (chr, str(peak[0]), str(peak[1]))
                    out.write(final)
                    peak_list.append(peak)


write()
