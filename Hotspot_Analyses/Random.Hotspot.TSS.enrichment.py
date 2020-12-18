#comparing hotspots and all genes (ensembl + isoseq annotations)
#replicating results of Hotspot paper (Shanfelter et al 2019) with new gene annotations
#this is the random enrichment script that will be run 10,000 times
#to run script: python3 Random.Hotspot.TSS.enrichment.py <transcript positions; combined Ensembl transcripts not in isoseq + Isoseq transcripts; converted to V5 positions> <random hotspot file with chr \t spot start \t spot end> <permutation number> <output file: perm number \t num hotspots at TSSs>
#Author: Alice Naftaly, Oct 2020

import sys

#read in transcript positions
#need to take strand into account
#returns dictionary with key == chr number and value = [transcript ID, TSS]
def read_transcripts():
    transcripts_file = sys.argv[1]
    transcript_dict = {}
    with open(transcripts_file ,'r') as transcripts:
        for line in transcripts:
            new_line = line.split()
            chr_num = new_line[0]
            transcript_id = new_line[3]
            strand = new_line[5]
            if strand == "+":
                tss = int(new_line[1])
            elif strand == "-":
                tss = int(new_line[2])
            dict_value = [transcript_id, tss]
            if chr_num in transcript_dict:
                transcript_dict[chr_num].append(dict_value)
            elif chr_num not in transcript_dict:
                transcript_dict.update({chr_num:[dict_value]})
    return transcript_dict

#read in random hotspots
#returns dictionary with key == chr number and value == [spot start, spot end]
def read_random_spots():
    spots_file = sys.argv[2]
    spots_dict = {}
    with open(spots_file, 'r') as spots:
        for line in spots:
            new_line = line.split()
            chr_num = new_line[0]
            spot_start = int(new_line[1])
            spot_end = int(new_line[2])
            dict_value = [spot_start, spot_end]
            if chr_num in spots_dict:
                spots_dict[chr_num].append(dict_value)
            elif chr_num not in spots_dict:
                spots_dict.update({chr_num:[dict_value]})
    return spots_dict

#count how many random hotspots fall within 3kb of TSS
#returns total number of hotspots that overlap TSS windows
def compare_TSS_spots():
    transcripts = read_transcripts()
    hotspots = read_random_spots()
    tss_overlap_total = []
    hotspot_overlapping_tss = []
    for chr in hotspots:
        single_chr_transcripts = transcripts[chr]
        single_chr_hotspots = hotspots[chr]
        for transcript in single_chr_transcripts:
            tss = transcript[1]
            upstream_tss = tss - 3000
            downstream_tss = tss + 3000
            transcript_id = transcript[0]
            for ind, hot in enumerate(single_chr_hotspots):
                hot_start = hot[0]
                hot_end = hot[1]
                if upstream_tss <= hot_start and upstream_tss <= hot_end and downstream_tss >= hot_start and downstream_tss >= hot_end:
                    hotspot_overlapping_tss.append(ind)
                elif upstream_tss > hot_start and upstream_tss < hot_end and downstream_tss > hot_start and downstream_tss > hot_end:
                    hotspot_overlapping_tss.append(ind)
                elif upstream_tss < hot_start and upstream_tss < hot_end and downstream_tss > hot_start and downstream_tss < hot_end:
                    hotspot_overlapping_tss.append(ind)
    set_hot_overlap = len(list(set(hotspot_overlapping_tss)))
    return set_hot_overlap

#write output
def write():
    overlap = compare_TSS_spots()
    perm_num = sys.argv[3]
    output = sys.argv[4]
    with open(output, 'a') as out:
        final = "%s\t%s\n" % (str(perm_num), str(overlap))
        out.write(final)

write()
