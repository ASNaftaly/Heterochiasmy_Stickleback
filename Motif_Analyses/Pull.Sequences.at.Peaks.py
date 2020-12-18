#pulling sequence under peaks
#position format:
#chr num \t peak start \t peak end
#need to read in one chromosome at a time, stickleback fasta, and peak positions
#will only examine autosomes
#to run script: python3 Pull.Sequences.at.Peaks.py <peak position file> <chr number> <stickleback fasta>  <output fasta file>
#Author: Alice Naftaly, Sept 2020


import sys

#read in peak positions file:
def read_peaks():
    peaks_file = sys.argv[1]
    chr = sys.argv[2]
    peaks_list = []
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                if chr == chr_num:
                    peak_start = int(new_line[1])
                    peak_end = int(new_line[2])
                    peak = [peak_start, peak_end]
                    peaks_list.append(peak)
    return peaks_list


#read stickleback fasta
def read_fasta():
    fasta_file = sys.argv[3]
    chr = sys.argv[2]
    fasta_dict = {}
    final_fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.split(" ")
                full_isoform_id = new_line[0].strip(" ")
                fasta_id = full_isoform_id.strip(">")
                final_fasta_id = fasta_id.strip("\n")
            else:
                new_line = line.strip("\n")
                fasta_dict.update({final_fasta_id:new_line})
    return fasta_dict[chr]


#pull sequences at peaks:
#returns dictionary with key == peak.startpos.endpos and value == sequence at peak
def pull_seqs():
    peaks = read_peaks()
    single_chr_fasta = read_fasta()
    peak_sequences = {}
    for peak in peaks:
        single_peak_header = "Peak.%s.%s" % (str(peak[0]), str(peak[1]))
        peak_start = peak[0]-1
        peak_end = peak[1]
        single_peak_sequence = single_chr_fasta[peak_start:peak_end]
        peak_sequences.update({single_peak_header:single_peak_sequence})
    return peak_sequences

#write output fasta
def write():
    peak_seqs = pull_seqs()
    chr_num = sys.argv[2]
    output = sys.argv[4]
    with open(output, 'a') as out:
        for peak in peak_seqs:
            header = ">" + chr_num + "." + peak + "\n"
            sequence = peak_seqs[peak]
            joined_seq = "".join(sequence)
            out.write(header)
            out.write(joined_seq + "\n")

write()
