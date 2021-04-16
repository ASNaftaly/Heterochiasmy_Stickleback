#Calculating the average read coverage per base pair for 1kb surrounding TSS across all chromosomes

#to run script: python3 Calc.Ave.TSS.Read.Cov.1kbwindow.py <all tss concatenated file> <output file>
#Author: Alice Naftaly, April 2021

import sys

#reading concatenated TSS file
#returns dictionary with key = chr number and value = read coverage per chromosome for all tsss per bp in 1kb around TSS
def read_all_tss():
    all_tss_file = sys.argv[1]
    read_cov_dict = {}
    with open(all_tss_file, 'r') as all_tss:
        for line in all_tss:
            new_line = line.split()
            chr_num = new_line[0]
            read_cov = int(new_line[2])
            if chr_num in read_cov_dict:
                read_cov_dict[chr_num].append(read_cov)
            elif chr_num not in read_cov_dict:
                read_cov_dict.update({chr_num:[read_cov]})
    return read_cov_dict


#calculate average read coverage per base pair
#returns dictionary with key = chr and value = average read coverage across all chromosomes
def calc_ave():
    read_cov_dict = read_all_tss()
    average_dict = {}
    for key in read_cov_dict:
        single_key = read_cov_dict[key]
        x = 0
        y = 1001
        while x < len(single_key):
            single_isoform = single_key[x:y]
            ave_single_isoform = round(sum(single_isoform)/1001, 2)
            if key in average_dict:
                average_dict[key].append(ave_single_isoform)
            elif key not in average_dict:
                average_dict.update({key:[ave_single_isoform]})
            x += 1001
            y += 1001
    return average_dict


#write average dictionary to file
def write():
    averages = calc_ave()
    output = sys.argv[2]
    with open(output, 'a') as out:
        chrs = ["chrI", "chrII","chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrXVII", "chrXVIII", "chrXX", "chrXXI"]
        for chr in chrs:
            single_key = averages[chr]
            for single in single_key:
                final = "%s\t%s\n" % (str(chr), str(single))
                out.write(final)

write()
