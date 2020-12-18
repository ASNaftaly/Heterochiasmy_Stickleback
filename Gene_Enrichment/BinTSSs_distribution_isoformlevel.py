#examine distribution of TSSs using Isoseq annotations
#will need to read in tss annotations as a bed file, make 10 bins per chromosome arm (will need centromere positions and chromosome sizes)

#to run script: python3 BinTSSs_distribution_isoformlevel.py <centromere positions file with roman numerals> <chr number as roman numeral> <chromosome sizes file> <TSS bed file> <output file>


import sys

#pulls centromere position from file based on chromosome number
#returns integer of centromere position in Mb
def pull_centromere():
    centromere_file = sys.argv[1]
    chr_num = sys.argv[2]
    with open(centromere_file, 'r') as centromeres:
        for line in centromeres:
            if line.startswith(chr_num + "\t"):
                centromere_pos = int(line.split()[1])
    return centromere_pos

#pulls length of chromosome from readcov file
#returns integer of chromosome length
def pull_chr_length():
    chr_length_file = sys.argv[3]
    chr_num = sys.argv[2]
    with open(chr_length_file, 'r') as readcov_file:
        for line in readcov_file:
            final_line = line.split()
            if chr_num == final_line[0]:
                chr_length = int(final_line[1])
    return chr_length

#create bins for chromosome
#arm 1 represents from the start of the chromosome to the centromere
#arm 2 represents from the centromere to the end of the chromosome
def create_bins():
    chr_length = pull_chr_length()
    centromere_pos = pull_centromere()
    arm_2 = chr_length - centromere_pos
    arm_1 = chr_length - arm_2
    bin_size_arm1 = int(round(arm_1/10,0))
    bin_size_arm2 = int(round(arm_2/10,0))
    #bins variable will have [bin size arm 1, bin size arm 2]
    bin_sizes = [bin_size_arm1,bin_size_arm2]
    return bin_sizes

#read TSS bed file
#bed file format = chr.num \t pos.one \t pos.two \t isoform id \t score \t strand
#need to pull strand, chr.num, and positions and isoform ids
#made sure to identify TSS based on strand
#returns dictionary with key== isoform id and value == [tss start, tss end]
def read_tss_positions():
    tss_bed_file = sys.argv[4]
    chr_num = sys.argv[2]
    tss_pos_dict = {}
    with open(tss_bed_file, 'r') as tss_bed:
        for line in tss_bed:
            new_line = line.split()
            if chr_num == new_line[0]:
                pos_one = new_line[1]
                pos_two = new_line[2]
                isoform_id = new_line[3]
                strand = new_line[5]
                if strand == "+":
                    tss_start = pos_one
                    tss_end = pos_two
                elif strand == "-":
                    tss_start = pos_two
                    tss_end = pos_one
                dict_value = [int(tss_start), int(tss_end)]
                tss_pos_dict.update({isoform_id:dict_value})
    return tss_pos_dict


#counting number of TSSs in bins for each arm
def count_tss_arm1():
    bins = create_bins()
    arm1_bin_size = bins[0]
    tss_pos = read_tss_positions()
    arm1_bin_dict = {}
    x = 0
    bin_start = 0
    bin_end = arm1_bin_size
    while x < 10:
        bin_count = 0
        for isoform in tss_pos:
            single_isoform = tss_pos[isoform]
            tss_start = single_isoform[0]
            tss_end = single_isoform[1]
            #if TSS is fully within the bin
            if tss_start >= bin_start and tss_start <= bin_end and tss_end >= bin_start and tss_end <= bin_end:
                bin_count += 1
            #if TSS starts in the bin but ends in the next bin
            elif tss_start >= bin_start and tss_start <= bin_end and tss_end >= bin_start and tss_end >= bin_end:
                bin_count += 1
            else:
                continue
        final = [bin_start,bin_end, bin_count]
        arm1_bin_dict.update({x:final})
        bin_start += arm1_bin_size
        bin_end += arm1_bin_size
        bin_count = 0
        x += 1
    return arm1_bin_dict

def count_tss_arm2():
    bins = create_bins()
    arm2_bin_size = bins[1]
    tss_pos = read_tss_positions()
    centromere_pos = pull_centromere()
    arm2_bin_dict = {}
    x = 10
    bin_start = centromere_pos
    bin_end = arm2_bin_size + centromere_pos
    while x < 20:
        bin_count = 0
        for isoform in tss_pos:
            single_isoform = tss_pos[isoform]
            tss_start = single_isoform[0]
            tss_end = single_isoform[1]
            #if TSS is fully within the bin
            if tss_start >= bin_start and tss_start <= bin_end and tss_end >= bin_start and tss_end <= bin_end:
                bin_count += 1
            #if TSS starts in the bin but ends in the next bin
            elif tss_start >= bin_start and tss_start <= bin_end and tss_end >= bin_start and tss_end >= bin_end:
                bin_count += 1
            else:
                continue
        final = [bin_start,bin_end, bin_count]
        arm2_bin_dict.update({x:final})
        bin_start += arm2_bin_size
        bin_end += arm2_bin_size
        bin_count = 0
        x += 1
    return arm2_bin_dict

#writing to one output file
def write():
    arm1_dict = count_tss_arm1()
    arm2_dict = count_tss_arm2()
    chr_num = sys.argv[2]
    bins = ["0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1.0"]
    output = sys.argv[5]
    with open(output,'a') as out:
        header = "Chr.Num\tBin.Num\tBin.Start\tBin.End\tNumber.Peaks\n"
        out.write(header)
        x = 0
        while x < 20:
            if x in arm1_dict:
                single_dict_value = arm1_dict[x]
                final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr_num),str(bins[x]), str(single_dict_value[0]),str(single_dict_value[1]),str(single_dict_value[2]))
                out.write(final)
                x += 1
            elif x in arm2_dict:
                single_dict_value = arm2_dict[x]
                final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr_num),str(bins[x]), str(single_dict_value[0]),str(single_dict_value[1]),str(single_dict_value[2]))
                out.write(final)
                x += 1

write()
