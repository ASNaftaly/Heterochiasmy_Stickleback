#pull 1kb around TSSs from autosomal TSS with readcoverage (4kb around TSS)
#will need to read in read coverage file (3 columns, chr num \t position (0-3999) \t read cov at bp)
#need to pull 1500 to 2500 for this window
#to run script: python3 Pull.1kb.window.TSSs.py <read coverage files with 4kb around TSS> <output file>
#Author: Alice Naftaly, March 2021

import sys

#read in read coverage file
#pulls just the 1kb around TSS
def read_file():
    rc_file = sys.argv[1]
    rc_dict = {}
    with open(rc_file, 'r') as rc:
        x = 0
        for line in rc:
            new_line = line.split()
            chr_num = new_line[0]
            position = int(new_line[1])
            read_cov = int(new_line[2])
            key = str(x) + "_" + chr_num
            if 1500 <= position <= 2500:
                if key in rc_dict:
                    rc_dict[key].append(read_cov)
                elif key not in rc_dict:
                    rc_dict.update({key:[read_cov]})
            elif position == 2501:
                x += 1
    return rc_dict

#calculating average read coverage per TSS
def calc_ave():
    read_cov = read_file()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for tss in read_cov:
            single_tss = read_cov[tss]
            ave_rc_single_tss = round(sum(single_tss)/len(single_tss),2)
            chr = tss.split("_")[1]
            final = "%s\t%s\n" % (str(chr), str(ave_rc_single_tss))
            out.write(final)


calc_ave()
