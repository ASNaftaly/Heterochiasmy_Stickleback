#need to pull isoforms from bed file and condense isoforms that have the same TSS so I can combine this with the average read coverage files
#to run script: python3 Condense.Isoforms.by.TSS.py < bed file with isoforms> <read coverage file for just autosomes> <output file>

import sys

#read in bed file:
def read_bed():
    bed_file = sys.argv[1]
    tss_dict = {}
    with open(bed_file, 'r') as bed:
        x = 0
        for line in bed:
            new_line = line.split()
            chr_num = new_line[0]
            tss_pos_1 = int(new_line[1])
            tss_pos_2 = int(new_line[2])
            isoform = new_line[3]
            strand = new_line[5]
            if strand == "+":
                tss = tss_pos_1
            elif strand == "-":
                tss = tss_pos_2
            dict_value = [isoform, tss, x]
            if chr_num in tss_dict:
                tss_dict[chr_num].append(dict_value)
            elif chr_num not in tss_dict:
                tss_dict.update({chr_num:[dict_value]})
            x += 1
    return tss_dict

#read average readcoverage file
def read_coverage():
    rc_file = sys.argv[2]
    rc_dict = {}
    with open(rc_file, 'r') as rc:
        for line in rc:
            new_line = line.split()
            chr_num = new_line[0]
            ave_read_cov = float(new_line[1])
            if chr_num in rc_dict:
                rc_dict[chr_num].append(ave_read_cov)
            elif chr_num not in rc_dict:
                rc_dict.update({chr_num:[ave_read_cov]})
    return rc_dict

#combine isoforms and average read coverage
def combine():
    rc_dict = read_coverage()
    isoforms = read_bed()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for chr in rc_dict:
            single_chr_rc = rc_dict[chr]
            single_chr_isoforms = isoforms[chr]
            tss_dict = {}
            x = 0
            for single_tss in single_chr_isoforms:
                iso = single_tss[0]
                tss_pos = single_tss[1]
                if tss_pos in tss_dict:
                    tss_dict[tss_pos].append(iso)
                elif tss_pos not in tss_dict:
                    tss_dict.update({tss_pos:[iso]})
            tss_list = []
            for val in tss_dict:
                single_val = tss_dict[val]
                tss_list.append(single_val)
            for index, v in enumerate(tss_list):
                read_cov = single_chr_rc[index]
                for a in v:
                    final = "%s\t%s\t%s\n" % (str(chr), str(a), str(read_cov))
                    out.write(final)


combine()
