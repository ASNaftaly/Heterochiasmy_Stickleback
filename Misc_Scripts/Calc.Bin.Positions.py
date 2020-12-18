#write bin positions to file
#to turn script: python3 Calc.Bin.Positions.py <centromere file> <chromosome sizes file>


import sys
#pulls centromere position from file based on chromosome number
def pull_centromere():
    centromere_file = sys.argv[1]
    centromere_positions = {}
    with open(centromere_file, 'r') as centromeres:
        for line in centromeres:
            if line.startswith("chr"):
                new_line = line.split()
                centromere_positions.update({new_line[0]:int(new_line[1])})
    return centromere_positions

#pulls length of chromosome from readcov file
#returns integer of chromosome length
def pull_chr_length():
    chr_length_file = sys.argv[2]
    chr_length = {}
    with open(chr_length_file, 'r') as readcov_file:
        for line in readcov_file:
            final_line = line.split()
            chr_length.update({final_line[0]:final_line[1]})
    return chr_length

#create bins for chromosome
#arm 1 represents from the start of the chromosome to the centromere
#arm 2 represents from the centromere to the end of the chromosome
def create_bins():
    chr_lengths = pull_chr_length()
    centromere_positions = pull_centromere()
    bins = {}
    for chr in chr_lengths:
        chr_length = int(chr_lengths[chr])
        centromere_pos = centromere_positions[chr]
        arm_2 = chr_length - centromere_pos
        arm_1 = chr_length - arm_2
        bin_size_arm1 = int(round(arm_1/10,0))
        bin_size_arm2 = int(round(arm_2/10,0))
        #bins variable will have [bin size arm 1, bin size arm 2]
        bin_sizes = [bin_size_arm1,bin_size_arm2]
        bins.update({chr:bin_sizes})
    return bins

#calc bin positions:
def calc_bins():
    bins = create_bins()
    centromeres = pull_centromere()
    bin_positions = {}
    for chr in bins:
        single_chr_bin = bins[chr]
        centromere = centromeres[chr]
        x = 0
        arm1_bin_size = single_chr_bin[0]
        bin_start = 0
        bin_end = bin_start + arm1_bin_size
        while x < 10:
            pos = [bin_start, bin_end]
            if chr in bin_positions:
                bin_positions[chr].append(pos)
            elif chr not in bin_positions:
                bin_positions.update({chr:[pos]})
            bin_start += arm1_bin_size
            bin_end += arm1_bin_size
            x += 1
        y = 10
        arm2_bin_size = single_chr_bin[1]
        bin_start_2 = centromere
        bin_end_2 = arm2_bin_size+bin_start_2
        while y < 20:
            pos = [bin_start_2, bin_end_2]
            if chr in bin_positions:
                bin_positions[chr].append(pos)
            elif chr not in bin_positions:
                bin_positions.update({chr:[pos]})
            bin_start_2 += arm2_bin_size
            bin_end_2 += arm2_bin_size
            y += 1
    return bin_positions

#write bin positions:
def write():
    bins = calc_bins()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for chr in bins:
            single_chr = bins[chr]
            for index, val in enumerate(single_chr):
                final = "%s\t%s\t%s\t%s\n" % (str(chr), str(index), str(val[0]), str(val[1]))
                out.write(final)


write()
