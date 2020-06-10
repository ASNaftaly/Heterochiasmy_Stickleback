#Pulling read coverage around TSSs using Isoseq annotations
#will do this for all TSSs
#this script is very similar to TSS_readcov_atac.py
#1. pull read coverage, transcript file in bed format
#2. assign upstream and downstream values for transcript = will start with 2kb from TSS (expect to see peak right at TSS that tapers off as distance from TSS increases)
#3. pull per bp read coverage for 4kb region around TSS (return dictionary with key = TSS and value = 4kb region with read coverage per bp)
#4. average read coverage per bp across all TSSs (i.e. combine all read coverages at position -2kb from TSS and take average)
#5. write average value per bp for all 4kb regions to a new file
#format:
#Position (0-4000) \t Ave.Read.Cov
#to run script: python3 Pull.ReadCov.all.TSSs.py <TSS positions in bed file> <chr # in roman numerals> <read coverage file with 1 chromosome> <output>
#author: Alice Naftaly, June 2020

import sys

#pulling TSSs from bed file
#returns list of transcription start sites
def pull_TSS():
    tss_file = sys.argv[1]
    chr_num = sys.argv[2]
    transcript_starts = []
    with open(tss_file,'r') as tss:
        for line in tss:
            new_line = line.split()
            chrnum = new_line[0]
            strand = new_line[5]
            if chrnum == chr_num:
                if strand == "+":
                    tss = int(new_line[1])
                    list_value = [strand, tss]
                    transcript_starts.append(list_value)
                elif strand == "-":
                    tss = int(new_line[2])
                    list_value = [strand, tss]
                    transcript_starts.append(list_value)
    return transcript_starts

#pull read coverage from file with single chromosome
#returns a dictionary with key = position and value = read coverage for that bp
def pull_read_cov():
    read_cov_file = sys.argv[3]
    read_cov_dict = {}
    with open(read_cov_file,'r') as read_coverage:
        for line in read_coverage:
            new_line = line.split()
            position = new_line[1]
            read_cov_per_bp = new_line[2]
            read_cov_dict.update({position:read_cov_per_bp})
    return read_cov_dict

#pull TSSs read coverage
#commented statements can be used to trouble shoot if there are any issues
#returns dictionary with key = tss position and value = [4001bp for plus/minus 2kb around TSS in order with [position, readcoverage] where position 2000 is the TSS]
def pull_tss_readcov():
    tss_positions = pull_TSS()
    read_coverage = pull_read_cov()
    tss_read_cov_dict = {}
    for tss in tss_positions:
        strand = tss[0]
        tss_pos = tss[1]
        tss_coverage = []
        if strand == "+":
            upstream_tss = tss_pos - 2000
            downstream_tss = tss_pos + 2000
            for position in read_coverage:
                if upstream_tss <= int(position) <= downstream_tss:
                    dict_value = [position, int(read_coverage[position])]
                    tss_coverage.append(dict_value)
                elif int(position) > downstream_tss:
                    break
            tss_read_cov_dict.update({tss_pos:tss_coverage})
        elif strand == "-":
            upstream_tss = tss_pos + 2000
            downstream_tss = tss_pos - 2000
            for position in read_coverage:
                if downstream_tss <= int(position) <= upstream_tss:
                    dict_value = [position, int(read_coverage[position])]
                    tss_coverage.append(dict_value)
                elif int(position) > upstream_tss:
                    break
            tss_coverage.reverse()
            tss_read_cov_dict.update({tss_pos:tss_coverage})
    return tss_read_cov_dict

#convert tss dict to be written
def convert_dict():
    tss_coverage = pull_tss_readcov()
    converted_tss_coverage = []
    for tss in tss_coverage:
        single_tss = tss_coverage[tss]
        x = 0
        while x < len(single_tss):
            single = single_tss[x]
            new_position = [x, single[1]]
            converted_tss_coverage.append(new_position)
            x += 1
    return converted_tss_coverage

#write output file
#writes coverage at every tss (4kb window) for each chromosome
def write():
    final_list = convert_dict()
    chr_num = sys.argv[2]
    output = sys.argv[4]
    with open(output, 'a') as out:
        for single_pos in final_list:
            final = "%s\t%s\t%s\n" % (str(chr_num),str(single_pos[0]), str(single_pos[1]))
            out.write(final)

write()
