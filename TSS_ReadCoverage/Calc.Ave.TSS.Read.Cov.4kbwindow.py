#Calculating the average read coverage per base pair for 4kb surrounding TSSs across all chromosomes
#script is very similar to Ave_tss_readcov_atac.py
#1. read TSS read coverage file for each sample (all chrs present except 19)
#2. calculate average read coverage per base pair
#3. write final 4kb to file
#to run script: python3 Calc.Ave.TSS.Read.Cov.4kbwindow.py <all tss concatenated file> <output file>
#Author: Alice Naftaly, June 2020

import sys

#reading concatenated TSS file
#returns dictionary with key = position (0-4000) and value = read coverage per chromosome for all tsss
def read_all_tss():
    all_tss_file = sys.argv[1]
    read_cov_dict = {}
    with open(all_tss_file, 'r') as all_tss:
        for line in all_tss:
            if line.startswith("Chr.Num"):
                continue
            else:
                new_line = line.split()
                position = new_line[1]
                read_cov = new_line[2]
                if position in read_cov_dict:
                    read_cov_dict[position].append(read_cov)
                elif position not in read_cov_dict:
                    read_cov_dict.update({position:[read_cov]})
    return read_cov_dict


#calculate average read coverage per base pair
#returns dictionary with key = position and value = average read coverage across all chromosomes
def calc_ave():
    read_cov_dict = read_all_tss()
    average_dict = {}
    for key in read_cov_dict:
        position = key
        single_position = read_cov_dict[position]
        int_single_position = [float(i) for i in single_position]
        ave_per_position = sum(int_single_position)/ len(int_single_position)
        average_dict.update({position:ave_per_position})
    return average_dict

#write average dictionary to file
def write():
    averages = calc_ave()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for key in averages:
            single_key = averages[key]
            final = "%s\t%s\n" % (str(key), str(single_key))
            out.write(final)

write()
