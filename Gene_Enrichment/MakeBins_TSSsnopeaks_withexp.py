#split peaks near TSSs into 20 bins across chromosomes to show if transcription has an effect
#1. read TSS file (peaks within 500bp of TSS), chromosome size file, centromere position file
#2. split chromosome into bins
#3. count the number of peaks near TSS in each bin
#4. write output to file
#to run script: python3 MakeBins_peaksnearTSS_atac.py <TSSs with no peaks but expressed file> <chr num> <centromere positions file> <Chr length file> <output file>
#Author: Alice Naftaly, August 2019, edited June 2020

import sys

#pull peaks at/near TSSs
#returns list of tss positions that are covered or surrounded by a peak
def read_tss_positions():
    tss_file = sys.argv[1]
    chr_num = sys.argv[2]
    tss_list = []
    iso_list = []
    with open(tss_file, 'r') as tss:
        for line in tss:
            if line.startswith(chr_num + "\t"):
                new_line = line.split()
                isoform = new_line[1]
                tss = int(new_line[2])
                if isoform in iso_list:
                    continue
                elif isoform not in iso_list:
                    iso_list.append(isoform)
                    tss_list.append(tss)
    return tss_list

#pulls centromere position from file based on chromosome number
#returns integer of centromere position in Mb
def pull_centromere():
    chr_num = sys.argv[2]
    centromere_file = sys.argv[3]
    centromere_pos = ""
    with open(centromere_file, 'r') as centromeres:
        for line in centromeres:
            if line.startswith(chr_num + "\t"):
                centromere_pos = int(line.split()[1])
    return centromere_pos

#pulls length of chromosome from readcov file
#returns integer of chromosome length
def pull_chr_length():
    chr_num = sys.argv[2]
    chr_length_file = sys.argv[4]
    chr_sizes_dict = {}
    with open(chr_length_file, 'r') as chr_sizes:
        for line in chr_sizes:
            new_line = line.split()
            chr_sizes_dict.update({new_line[0]:new_line[1]})
    return chr_sizes_dict[chr_num]


#create bins for chromosome
#arm 1 represents from the start of the chromosome to the centromere
#arm 2 represents from the centromere to the end of the chromosome
def create_bins():
    chr_length = pull_chr_length()
    centromere_pos = pull_centromere()
    arm_2 = int(chr_length) - centromere_pos
    arm_1 = int(chr_length) - arm_2
    bin_size_arm1 = int(round(arm_1/10,0))
    bin_size_arm2 = int(round(arm_2/10,0))
    #bins variable will have [bin size arm 1, bin size arm 2]
    bin_sizes = [bin_size_arm1,bin_size_arm2]
    return bin_sizes


#count tss in arm 1
def count_tss_arm1():
    bins = create_bins()
    arm1_bin_size = bins[0]
    tss = read_tss_positions()
    bin_dict = {}
    x = 0
    bin_start = 0
    bin_end = arm1_bin_size
    while x < 10:
        bin_count = 0
        for index, value in enumerate(tss):
            if bin_start <= int(value) <= bin_end:
                bin_count += 1
            else:
                continue
        final = [bin_start, bin_end, bin_count]
        bin_dict.update({x:final})
        bin_start += arm1_bin_size
        bin_end += arm1_bin_size
        bin_count = 0
        x += 1
    return bin_dict

#count tss in arm 2
def count_tss_arm2():
    bins = create_bins()
    centromere = pull_centromere()
    chr_size = int(pull_chr_length())
    arm2_bin_size = bins[1]
    tss = read_tss_positions()
    bin_dict = {}
    x = 10
    bin_start = centromere
    bin_end = arm2_bin_size + centromere
    while x < 19:
        bin_count = 0
        for index, value in enumerate(tss):
            if bin_start <= int(value) <= bin_end:
                bin_count += 1
            else:
                continue
        final = [bin_start, bin_end, bin_count]
        bin_dict.update({x:final})
        bin_start += arm2_bin_size
        bin_end += arm2_bin_size
        bin_count = 0
        x += 1
    for index, value in enumerate(tss):
        if bin_start <= int(value):
            bin_count += 1
        else:
            continue
    final = [bin_start, chr_size, bin_count]
    bin_dict.update({19:final})
    return bin_dict

#writing to one output file
def write():
    arm1_dict = count_tss_arm1()
    arm2_dict = count_tss_arm2()
    bins = ["0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1.0"]
    output = sys.argv[5]
    with open(output,'a') as out:
        header = "Bin.Num\tBin.Start\tBin.End\tNumber.TSSs\n"
        out.write(header)
        x = 0
        while x < 20:
            if x in arm1_dict:
                single_dict_value = arm1_dict[x]
                final = "%s\t%s\t%s\t%s\n" % (str(bins[x]), str(single_dict_value[0]),str(single_dict_value[1]),str(single_dict_value[2]))
                out.write(final)
                x += 1
            elif x in arm2_dict:
                single_dict_value = arm2_dict[x]
                final = "%s\t%s\t%s\t%s\n" % (str(bins[x]), str(single_dict_value[0]),str(single_dict_value[1]),str(single_dict_value[2]))
                out.write(final)
                x += 1

write()
