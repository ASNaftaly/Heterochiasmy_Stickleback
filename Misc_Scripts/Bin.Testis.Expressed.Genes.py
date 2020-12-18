#bin isoforms that are testis biased or testis expressed into bins
#to run script: python3 Bin.Testis.Expressed.Genes.py


import sys

#read in isoforms:
def read_isoforms():
    isoform_file = sys.argv[1]
    all_isoforms = []
    with open(isoform_file ,'r') as isoforms:
        for line in isoforms:
            new_line = line.strip("\n")
            all_isoforms.append(new_line)
    return all_isoforms


#read bed file and pull start position
def pull_start_positions():
    isoforms = read_isoforms()
    bed_file = sys.argv[2]
    start_pos_dict = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            new_line = line.split()
            chr = new_line[0]
            isoform_id = new_line[3]
            strand = new_line[5]
            if strand == "+":
                start = int(new_line[1])
            elif start == "-":
                start = int(new_line[2])
            if isoform_id in isoforms:
                if chr in start_pos_dict:
                    start_pos_dict[chr].append(start)
                elif chr not in start_pos_dict:
                    start_pos_dict.update({chr:[start]})
    return start_pos_dict

#read bin_positions
def read_bin_positions():
    bin_file = sys.argv[3]
    bin_dict = {}
    with open(bin_file, 'r') as bins:
        for line in bins:
            new_line = line.split()
            chr = new_line[0]
            bin_num = new_line[1]
            bin_start = int(new_line[2])
            bin_end = int(new_line[3])
            dict_value = [bin_num, bin_start, bin_end]
            if chr in bin_dict:
                bin_dict[chr].append(dict_value)
            elif chr not in bin_dict:
                bin_dict.update({chr:[dict_value]})
    return bin_dict


#bin isoforms
def bin_isoforms():
    isoform_starts = pull_start_positions()
    bin_positions = read_bin_positions()
    binned_counts = {}
    for chr in isoform_starts:
        single_chr_isoforms = isoform_starts[chr]
        single_chr_bins = bin_positions[chr]
        for bin in single_chr_bins:
            bin_num = bin[0]
            bin_start = bin[1]
            bin_end = bin[2]
            bin_count = 0
            for start in single_chr_isoforms:
                if bin_start <= start <= bin_end:
                    bin_count += 1
            if bin_num in binned_counts:
                binned_counts[bin_num].append(bin_count)
            elif bin_num not in binned_counts:
                binned_counts.update({bin_num:[bin_count]})
    return binned_counts

#total bins:
def 
