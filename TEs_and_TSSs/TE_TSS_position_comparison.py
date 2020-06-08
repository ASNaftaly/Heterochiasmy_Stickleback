#comparing locations of annotated TEs with annotated TSS (using Isoseq annotations)
#this script will be similar to TEs_near_TSS.py
#will need to read in the bed file with the correct TSS to HiC assembly (or whatever current assembly is being used) and the repeat masker GTF file
#will be excluding the X and Y chromosomes (chrXIX) because the annotations on this chromosome might be a little off given the similarity between the X and Y
#to run script: python3 TE_TSS_position_comparison.py <bed file with all transcript positions> <repeat masker gtf file> <output file; summary file> <output file, tss overlap> <output file extended tss overlap> <output file, tts overlap> <output file, extended tts overlap>
#author: Alice Naftaly, June 2020

import sys

#read in bed file
#returns dictionary with key == chr number and value == [isoform id, 5' position on + strand, 3' position on + strand, strand]
def read_isoform_bed():
    bed_file = sys.argv[1]
    bed_dict = {}
    with open(bed_file, 'r') as bed_info:
        for line in bed_info:
            new_line = line.split("\t")
            chr_num = new_line[0]
            start_pos = int(new_line[1])
            end_pos = int(new_line[2])
            isoform_id = new_line[3]
            strand = new_line[5].strip("\n")
            dict_value = [isoform_id, start_pos, end_pos, strand]
            if chr_num in bed_dict:
                bed_dict[chr_num].append(dict_value)
            elif chr_num not in bed_dict:
                bed_dict.update({chr_num:[dict_value]})
    return bed_dict


#read in repeat masker GTF
#returns dictionary with key == DNA, LINE, SINE, SINE?, LTR, RC, low complexity, satellite, Simple repeat, snRNA, or Unknown and value == [repeat name, start pos, end pos, strand]
def read_repeat_gtf():
    gtf_file = sys.argv[2]
    repeat_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            chr_num = new_line[0]
            broad_te_class = new_line[2]
            start_pos = int(new_line[3])
            end_pos = int(new_line[4])
            strand = new_line[6]
            repeat_info = new_line[8].split(" ")
            repeat_name_full = repeat_info[3].strip(";")
            repeat_name_final = repeat_name_full.strip("\'")
            dict_value = [repeat_name_final, start_pos, end_pos, strand, chr_num]
            if broad_te_class in repeat_dict:
                repeat_dict[broad_te_class].append(dict_value)
            elif broad_te_class not in repeat_dict:
                repeat_dict.update({broad_te_class:[dict_value]})
    return repeat_dict

#for comparing positions of TEs vs TSS with do this in 2 ways:
#1: Does TE overlap TSS (without adjusting the TSS position at all)
#   on same strand or different strand?
#2: Does TE overlap TSS +/- 500bp
#   on same strand or different strand
def count_TSS_overlap():
    tss_positions = read_isoform_bed()
    te_positions = read_repeat_gtf()
    tss_te_overlap = {}
    for repeat_type in te_positions:
        single_repeat_type = te_positions[repeat_type]
        for single_repeat in single_repeat_type:
            repeat_name = single_repeat[0]
            repeat_start = single_repeat[1]
            repeat_end = single_repeat[2]
            repeat_strand = single_repeat[3]
            repeat_chr_num = single_repeat[4]
            if repeat_chr_num in tss_positions:
                single_chr_positions = tss_positions[repeat_chr_num]
                for single_tss in single_chr_positions:
                    isoform_id = single_tss[0]
                    isoform_start = single_tss[1]
                    isoform_strand = single_tss[3]
                    if repeat_start <= isoform_start <= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
    return tss_te_overlap

def count_TSS_extended_overlap():
    tss_positions = read_isoform_bed()
    te_positions = read_repeat_gtf()
    tss_te_overlap = {}
    for repeat_type in te_positions:
        single_repeat_type = te_positions[repeat_type]
        for single_repeat in single_repeat_type:
            repeat_name = single_repeat[0]
            repeat_start = single_repeat[1]
            repeat_end = single_repeat[2]
            repeat_strand = single_repeat[3]
            repeat_chr_num = single_repeat[4]
            if repeat_chr_num in tss_positions:
                single_chr_positions = tss_positions[repeat_chr_num]
                for single_tss in single_chr_positions:
                    isoform_id = single_tss[0]
                    isoform_start = single_tss[1]
                    isoform_strand = single_tss[3]
                    upstream_tss = isoform_start - 500
                    downstream_tss = isoform_start + 500
                    if repeat_start <= isoform_start <= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                    #if TE is within 1kb window
                    elif upstream_tss <= repeat_start and repeat_start <= downstream_tss and repeat_end >= upstream_tss and downstream_tss >= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                    #if TE overlaps upstream portion of TSS window:
                    elif repeat_start <= upstream_tss and repeat_start <= downstream_tss and repeat_end >= upstream_tss and downstream_tss >= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                    #if TE overlaps downstream_portion of TSS window:
                    elif repeat_start >= upstream_tss and repeat_start <= downstream_tss and repeat_end >= upstream_tss and downstream_tss <= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tss_te_overlap:
                                tss_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tss_te_overlap:
                                tss_te_overlap.update({repeat_type:[value]})
    return tss_te_overlap


def count_TTS_overlap():
    tss_positions = read_isoform_bed()
    te_positions = read_repeat_gtf()
    tts_te_overlap = {}
    for repeat_type in te_positions:
        single_repeat_type = te_positions[repeat_type]
        for single_repeat in single_repeat_type:
            repeat_name = single_repeat[0]
            repeat_start = single_repeat[1]
            repeat_end = single_repeat[2]
            repeat_strand = single_repeat[3]
            repeat_chr_num = single_repeat[4]
            if repeat_chr_num in tss_positions:
                single_chr_positions = tss_positions[repeat_chr_num]
                for single_tss in single_chr_positions:
                    isoform_id = single_tss[0]
                    isoform_end = single_tss[2]
                    isoform_strand = single_tss[3]
                    if repeat_start <= isoform_end <= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
    return tts_te_overlap

def count_TTS_extended_overlap():
    tss_positions = read_isoform_bed()
    te_positions = read_repeat_gtf()
    tts_te_overlap = {}
    for repeat_type in te_positions:
        single_repeat_type = te_positions[repeat_type]
        extended_overlap_same_strand = []
        extended_overlap_diff_strand = []
        for single_repeat in single_repeat_type:
            repeat_name = single_repeat[0]
            repeat_start = single_repeat[1]
            repeat_end = single_repeat[2]
            repeat_strand = single_repeat[3]
            repeat_chr_num = single_repeat[4]
            if repeat_chr_num in tss_positions:
                single_chr_positions = tss_positions[repeat_chr_num]
                for single_tss in single_chr_positions:
                    isoform_id = single_tss[0]
                    isoform_end = single_tss[2]
                    isoform_strand = single_tss[3]
                    upstream_tts = isoform_end - 500
                    downstream_tts = isoform_end + 500
                    if repeat_start <= isoform_end <= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                    #if TE is within 1kb window
                    elif upstream_tts <= repeat_start and repeat_start <= downstream_tts and repeat_end >= upstream_tts and downstream_tts >= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                    #if TE overlaps upstream portion of TTS window:
                    elif repeat_start <= upstream_tts and repeat_start <= downstream_tts and repeat_end >= upstream_tts and downstream_tts >= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                    #if TE overlaps downstream_portion of TTS window:
                    elif repeat_start >= upstream_tts and repeat_start <= downstream_tts and repeat_end >= upstream_tts and downstream_tts <= repeat_end:
                        if isoform_strand == repeat_strand:
                            value = ["same_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
                        elif isoform_strand != repeat_strand:
                            value = ["diff_strand", isoform_id, repeat_name]
                            if repeat_type in tts_te_overlap:
                                tts_te_overlap[repeat_type].append(value)
                            elif repeat_type not in tts_te_overlap:
                                tts_te_overlap.update({repeat_type:[value]})
    return tts_te_overlap


#create summary counts separately for each overlap
def create_summary_tss_overlap():
    tss_overlap = count_TSS_overlap()
    summary_dict = {}
    TE_types = ["DNA", "LINE", "SINE", "SINE?", "LTR", "RC", "Low_complexity", "Satellite", "Simple_repeat", "snRNA", "Unknown"]
    for te in TE_types:
        if te in tss_overlap:
            single_te_type = tss_overlap[te]
            total_counts = len(single_te_type)
            same_strand_counts = 0
            diff_strand_counts = 0
            for single in single_te_type:
                if single[0] == "same_strand":
                    same_strand_counts += 1
                elif single[0] == "diff_strand":
                    diff_strand_counts += 1
            final_summary = [str(total_counts), str(same_strand_counts), str(diff_strand_counts)]
            summary_dict.update({te:final_summary})
        elif te not in tss_overlap:
            final_summary = ["0", "0", "0"]
            summary_dict.update({te:final_summary})
    return summary_dict

def create_summary_extended_tss_overlap():
    tss_overlap = count_TSS_extended_overlap()
    summary_dict = {}
    TE_types = ["DNA", "LINE", "SINE", "SINE?", "LTR", "RC", "Low_complexity", "Satellite", "Simple_repeat", "snRNA", "Unknown"]
    for te in TE_types:
        if te in tss_overlap:
            single_te_type = tss_overlap[te]
            total_counts = len(single_te_type)
            same_strand_counts = 0
            diff_strand_counts = 0
            for single in single_te_type:
                if single[0] == "same_strand":
                    same_strand_counts += 1
                elif single[0] == "diff_strand":
                    diff_strand_counts += 1
            final_summary = [str(total_counts), str(same_strand_counts), str(diff_strand_counts)]
            summary_dict.update({te:final_summary})
        elif te not in tss_overlap:
            final_summary = ["0", "0", "0"]
            summary_dict.update({te:final_summary})
    return summary_dict

def create_summary_tts_overlap():
    tts_overlap = count_TTS_overlap()
    summary_dict = {}
    TE_types = ["DNA", "LINE", "SINE", "SINE?", "LTR", "RC", "Low_complexity", "Satellite", "Simple_repeat", "snRNA", "Unknown"]
    for te in TE_types:
        if te in tts_overlap:
            single_te_type = tts_overlap[te]
            total_counts = len(single_te_type)
            same_strand_counts = 0
            diff_strand_counts = 0
            for single in single_te_type:
                if single[0] == "same_strand":
                    same_strand_counts += 1
                elif single[0] == "diff_strand":
                    diff_strand_counts += 1
            final_summary = [str(total_counts), str(same_strand_counts), str(diff_strand_counts)]
            summary_dict.update({te:final_summary})
        elif te not in tts_overlap:
            final_summary = ["0", "0", "0"]
            summary_dict.update({te:final_summary})
    return summary_dict

def create_summary_tts_overlap():
    tts_overlap = count_TTS_extended_overlap()
    summary_dict = {}
    TE_types = ["DNA", "LINE", "SINE", "SINE?", "LTR", "RC", "Low_complexity", "Satellite", "Simple_repeat", "snRNA", "Unknown"]
    for te in TE_types:
        if te in tts_overlap:
            single_te_type = tts_overlap[te]
            total_counts = len(single_te_type)
            same_strand_counts = 0
            diff_strand_counts = 0
            for single in single_te_type:
                if single[0] == "same_strand":
                    same_strand_counts += 1
                elif single[0] == "diff_strand":
                    diff_strand_counts += 1
            final_summary = [str(total_counts), str(same_strand_counts), str(diff_strand_counts)]
            summary_dict.update({te:final_summary})
        elif te not in tts_overlap:
            final_summary = ["0", "0", "0"]
            summary_dict.update({te:final_summary})
    return summary_dict


#write summary file
def write_summary():
    tss_overlap = create_summary_tss_overlap()
    extended_tss_overlap = create_summary_extended_tss_overlap()
    tts_overlap = create_summary_tts_overlap()
    extended_tts_overlap = create_summary_extended_tss_overlap()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "TE.Type\tTotal.TSS.Overlap\tTSS.Overlap.Same.Strand\tTSS.Overlap.Diff.Strand\tTotal.Extended.TSS.Overlap\tExtended.TSS.Overlap.Same.Strand\tExtended.TSS.Overlap.Diff.Strand\tTotal.TTS.Overlap\tTTS.Overlap.Same.Strand\tTTS.Overlap.Diff.Strand\tTotal.Extended.TTS.Overlap\tExtended.TTS.Overlap.Same.Strand\tExtended.TTS.Overlap.Diff.Strand\n"
        out.write(header)
        for te in tss_overlap:
            single_tss_overlap = tss_overlap[te]
            single_extended_tss_overlap = extended_tss_overlap[te]
            single_tts_overlap = tts_overlap[te]
            single_extended_tts_overlap = extended_tts_overlap[te]
            final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(te), str(single_tss_overlap[0]), str(single_tss_overlap[1]), str(single_tss_overlap[2]), str(single_extended_tss_overlap[0]), str(single_extended_tss_overlap[1]), str(single_extended_tss_overlap[2]), str(single_tts_overlap[0]), str(single_tts_overlap[1]), str(single_tts_overlap[2]), str(single_extended_tts_overlap[0]), str(single_extended_tts_overlap[1]), str(single_extended_tts_overlap[2]))
            out.write(final)

#write output for each overlap set
def write_tss_overlap():
    tss_overlap = count_TSS_overlap()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Broad.TE.Class\tIsoform.ID\tRepeat.Name\tSame.or.Diff.Strand\n"
        out.write(header)
        for te in tss_overlap:
            single_te = tss_overlap[te]
            for single in single_te:
                strand = single[0]
                isoform = single[1]
                repeat = single[2]
                if strand == "same_strand":
                    final_strand = "same"
                elif strand == "diff_strand":
                    final_strand = "diff"
                final = "%s\t%s\t%s\t%s\n" % (str(te), str(isoform), str(repeat), str(final_strand))
                out.write(final)

def write_extended_tss_overlap():
    tss_overlap = count_TSS_extended_overlap()
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = "Broad.TE.Class\tIsoform.ID\tRepeat.Name\tSame.or.Diff.Strand\n"
        out.write(header)
        for te in tss_overlap:
            single_te = tss_overlap[te]
            for single in single_te:
                strand = single[0]
                isoform = single[1]
                repeat = single[2]
                if strand == "same_strand":
                    final_strand = "same"
                elif strand == "diff_strand":
                    final_strand = "diff"
                final = "%s\t%s\t%s\t%s\n" % (str(te), str(isoform), str(repeat), str(final_strand))
                out.write(final)

def write_tts_overlap():
    tts_overlap = count_TTS_overlap()
    output = sys.argv[6]
    with open(output, 'a') as out:
        header = "Broad.TE.Class\tIsoform.ID\tRepeat.Name\tSame.or.Diff.Strand\n"
        out.write(header)
        for te in tts_overlap:
            single_te = tts_overlap[te]
            for single in single_te:
                strand = single[0]
                isoform = single[1]
                repeat = single[2]
                if strand == "same_strand":
                    final_strand = "same"
                elif strand == "diff_strand":
                    final_strand = "diff"
                final = "%s\t%s\t%s\t%s\n" % (str(te), str(isoform), str(repeat), str(final_strand))
                out.write(final)

def write_extended_tts_overlap():
    tts_overlap = count_TTS_extended_overlap()
    output = sys.argv[7]
    with open(output, 'a') as out:
        header = "Broad.TE.Class\tIsoform.ID\tRepeat.Name\tSame.or.Diff.Strand\n"
        out.write(header)
        for te in tts_overlap:
            single_te = tts_overlap[te]
            for single in single_te:
                strand = single[0]
                isoform = single[1]
                repeat = single[2]
                if strand == "same_strand":
                    final_strand = "same"
                elif strand == "diff_strand":
                    final_strand = "diff"
                final = "%s\t%s\t%s\t%s\n" % (str(te), str(isoform), str(repeat), str(final_strand))
                out.write(final)


#call all functions()
def call():
    summary = write_summary()
    tss_overlap = write_tss_overlap()
    extended_tss_overlap = write_extended_tss_overlap()
    tts_overlap = write_tts_overlap()
    extended_tts_overlap = write_extended_tts_overlap()

call()
