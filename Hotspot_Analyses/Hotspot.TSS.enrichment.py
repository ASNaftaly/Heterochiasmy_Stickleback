#comparing hotspots and all genes (ensembl + isoseq annotations)
#replicating results of Hotspot paper (Shanfelter et al 2019) with new gene annotations
#will ask how many hotspots fall within 3kb of a TSS
#to run script: python3 Hotspot.TSS.enrichment.py <transcript positions; combined Ensembl transcripts not in isoseq + Isoseq transcripts; converted to V5 positions> <hotspots file in bed format; converted to V5 positions> <output file; has one hotspot per line that is within 3kb of a hotspot> >> <standard output has total number of TSSs with hotspots in window>
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

#read in hotspots
#returns dictionary with key == chr number and value == [hotspot start, hotspot end, hotspots ID]
def read_hotspots():
    hotspots_file = sys.argv[2]
    hotspots_dict = {}
    with open(hotspots_file, 'r') as hotspots:
        for line in hotspots:
            new_line = line.split()
            chr_num = new_line[0]
            hot_start = int(new_line[1])
            hot_end = int(new_line[2])
            hot_id = new_line[3]
            dict_value = [hot_start, hot_end, hot_id]
            if chr_num in hotspots_dict:
                hotspots_dict[chr_num].append(dict_value)
            elif chr_num not in hotspots_dict:
                hotspots_dict.update({chr_num:[dict_value]})
    return hotspots_dict

#count how many hotspots fall within 3kb of TSS
#returns total number of TSSs within 3kb of TSSs in a list, and list of hotspots that overlap TSS windows
def compare_TSS_hotspots():
    transcripts = read_transcripts()
    hotspots = read_hotspots()
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
            tss_overlap = 0
            for hot in single_chr_hotspots:
                hot_start = hot[0]
                hot_end = hot[1]
                hot_id = hot[2]
                if upstream_tss <= hot_start and upstream_tss <= hot_end and downstream_tss >= hot_start and downstream_tss >= hot_end:
                    tss_overlap += 1
                    hotspot_overlapping_tss.append(hot_id)
                elif upstream_tss > hot_start and upstream_tss < hot_end and downstream_tss > hot_start and downstream_tss > hot_end:
                    tss_overlap += 1
                    hotspot_overlapping_tss.append(hot_id)
                elif upstream_tss < hot_start and upstream_tss < hot_end and downstream_tss > hot_start and downstream_tss < hot_end:
                    tss_overlap += 1
                    hotspot_overlapping_tss.append(hot_id)
            if tss_overlap > 0:
                tss_overlap_total.append(tss_overlap)
    sum_tss_overlap = sum(tss_overlap_total)
    set_hot_overlap = list(set(hotspot_overlapping_tss))
    return sum_tss_overlap, set_hot_overlap

#write output
def write():
    tss_overlap, hot_overlap = compare_TSS_hotspots()
    final_tss_overlap = "Total TSSs that have hotspots within 3kb: %s\n" % str(tss_overlap)
    print(final_tss_overlap)
    output = sys.argv[3]
    with open(output, 'a') as out:
        for hot in hot_overlap:
            final_hot = "%s\n" % str(hot)
            out.write(final_hot)

write()
