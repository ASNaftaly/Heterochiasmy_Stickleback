#filter out peaks that do not meet criteria specified by manual filters
##will need to read in average read depth from bed graph file and xls peak files
#to run script: python3 FilterPeaks.py <bedgraph file> <xls file> <read coverage minimum> <pileup minimum> <output file>
#Author: Alice Naftaly, Jan 2022


import sys


#read in bedgraph read coverage file
#just need chr, peak start, and normalized read counts
#returns dictionary with key == chr number and value == [peak start position, normalized read counts]
def read_bedgraph():
    bedgraph_file = sys.argv[1]
    bedgraph_dict = {}
    with open(bedgraph_file, 'r') as bedgraph:
        for line in bedgraph:
            new_line = line.split("\t")
            chr_num = new_line[0]
            peak_start = int(new_line[1])
            read_count = int(new_line[2])
            dict_value = [peak_start, read_count]
            if str(chr_num) in bedgraph_dict:
                bedgraph_dict[str(chr_num)].append(dict_value)
            elif str(chr_num) not in bedgraph_dict:
                bedgraph_dict.update({str(chr_num):[dict_value]})
    return bedgraph_dict

#read peak file
#need whole line but will need access to peak start (column 2) and pileup (column 6)
#returns dictionary with key == chr number and value == full line split on tabs
def read_peaks():
    peaks_file = sys.argv[2]
    peaks_dict = {}
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("chr"):
                new_line = line.split()
                if new_line[1] == "start":
                    continue
                else:
                    chr_num = new_line[0]
                    if chr_num in peaks_dict:
                        peaks_dict[chr_num].append(new_line)
                    elif chr_num not in peaks_dict:
                        peaks_dict.update({chr_num:[new_line]})
    return peaks_dict

#remove peaks that do not meet filter
#returns dictionary with key == chr and value == line from xls file
def filter_peaks():
    peaks = read_peaks()
    normalized_read_coverage = read_bedgraph()
    min_norm_read_cov = int(sys.argv[3])
    min_pileup = int(sys.argv[4])
    filtered_peaks = {}
    for chr in peaks:
        single_chr_peaks = peaks[chr]
        single_chr_norm = normalized_read_coverage[chr]
        x = 0
        while x < len(single_chr_peaks):
            norm_coverage = single_chr_norm[x][1]
            pileup = float(single_chr_peaks[x][5])
            if norm_coverage > min_norm_read_cov and pileup > min_pileup:
                if chr in filtered_peaks:
                    filtered_peaks[chr].append(single_chr_peaks[x])
                elif chr not in filtered_peaks:
                    filtered_peaks.update({chr:[single_chr_peaks[x]]})
            x += 1
    return filtered_peaks


#write new xls output file
def write():
    filtered_peaks = filter_peaks()
    output = sys.argv[5]
    chrs = ["chrI","chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrXVII", "chrXVIII", "chrXX", "chrXXI"]
    with open(output, 'a') as out:
        for chr in chrs:
            single_chr_peaks = filtered_peaks[chr]
            for single in single_chr_peaks:
                final = "\t".join(single)
                out.write(final + "\n")



write()
