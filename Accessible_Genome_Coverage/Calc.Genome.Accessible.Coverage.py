#calculating the amount of the genome covered by peaks
#can do this with any peak file in the xls format
#this will not take into account peak number and will not be normalized (need to figure out the most appropriate way to normalize based on number of reads)
#output format = Chr.Num \t Chr.Size \t Bp.Covered.by.ATAC.peak \t Percent.Covered.by.ATAC.peak \n
#to run script: python3 Calc.Genome.Accessible.Coverage.py <peak file xls> <chr sizes file> <output file>
#Author: Alice Naftaly, June 2020


import sys

#read in peak file
#returns dictionary with key == chr number and value == list of peak sizes
def read_peaks():
    peaks_file = sys.argv[1]
    peaks_dict = {}
    with open(peaks_file, 'r') as peaks:
        for line in peaks:
            if line.startswith("chr\tstart"):
                continue
            elif line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                peak_size = int(new_line[3])
                if chr_num in peaks_dict:
                    peaks_dict[chr_num].append(peak_size)
                elif chr_num not in peaks_dict:
                    peaks_dict.update({chr_num:[peak_size]})
    return peaks_dict

#read in chr sizes
#returns dictionary with key == chr num and value == chr size
def read_chr_sizes():
    chr_size_file = sys.argv[2]
    chr_size_dict = {}
    with open(chr_size_file, 'r') as sizes:
        for line in sizes:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                chr_size = int(new_line[1])
                chr_size_dict.update({chr_num:chr_size})
    return chr_size_dict


#calculating total chr covered by ATAC peaks
#returns dictionary with key == chr num and value == bps covered by ATAC peaks
def calc_coverage():
    peaks_dict = read_peaks()
    coverage_dict = {}
    for chr in peaks_dict:
        single_chr = peaks_dict[chr]
        sum_single_chr = sum(single_chr)
        coverage_dict.update({chr:sum_single_chr})
    return coverage_dict

#write un-normalized summary file
def write():
    chr_sizes = read_chr_sizes()
    peak_coverage = calc_coverage()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Chr.Num\tChr.Size\tBp.Covered.by.ATAC.peak\tPercent.Covered.by.ATAC.peaks\n"
        out.write(header)
        for chr in chr_sizes:
            if chr in peak_coverage:
                chr_size = chr_sizes[chr]
                single_chr_coverage = peak_coverage[chr]
                percent_coverage = round((single_chr_coverage/chr_size)*100,2)
                final = "%s\t%s\t%s\t%s\n" % (str(chr), str(chr_size), str(single_chr_coverage), str(percent_coverage))
                out.write(final)


write()
