#Pulling Transcript IDS for genes enriched around Peaks
#Peaks reported if within 500bp of a TSS
#will complete for both tissues and samples
#to run script: python3 PullGeneIDs_atac.py <TSS observed file> <chr number in roman numerals> <TSS file> <output file>
#Author: Alice Naftaly, July 2019, edited August 2019

'''to run script:
see PullGeneIDs_*.sh
'''

import sys

#pull TSS position from enrichment files
#returns a list of TSS positions
def pull_peaks_at_TSS():
    tss_observed_file = sys.argv[1]
    peak_tss_pos = []
    with open(tss_observed_file, 'r') as tss_obs:
        for line in tss_obs:
            new_line = line.split()
            if new_line[0] == "Peak.Start":
                continue
            else:
                peak_tss_pos.append(new_line[1])
    return peak_tss_pos


#pulling transcripts from TSS file
#need chromosome number
#returns as a dictionary
#key = tss start site
#value = tss ID for ensembl
def Pull_TSS_info():
    data = []
    chr_num = sys.argv[2]
    psl_file = sys.argv[3]
    tss_dict = {}
    with open (psl_file,'r') as glaz_file:
        for line in glaz_file:
            new_line = line.split()
            if new_line[13] == chr_num:
                tss_start = new_line[15]
                tss_ID = new_line[9]
                tss_dict.update({tss_start:tss_ID})
    return tss_dict


#Comparing TSSs enriched at hotspots to TSS info list for each chromosome
#Will pull tss and ID
#This does not create duplicates (if there is more than 1 peak near TSS, it is only counted once)
#final file has Chr.Num \t TSS.pos \t Transcript.ID
def pull_ID():
    peak_at_tss = pull_peaks_at_TSS()
    tss_info = Pull_TSS_info()
    chr_num = sys.argv[2]
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Chr.Num\tTSS\tTranscrip.ID\n"
        out.write(header)
        x = 0
        for key in tss_info:
            if key in peak_at_tss:
                tss = key
                gene_id = tss_info[key]
                final = "%s\t%s\t%s\n" % (str(chr_num),str(tss),str(gene_id))
                out.write(final)

pull_ID()
