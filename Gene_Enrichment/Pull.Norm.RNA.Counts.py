#script to pull normalized RNA read counts for TSSs that have peaks
#will read in observed TSS overlap file (has tss information and peak position), normalized counts csv, isoform to gene id conversion file
#output file format:
#chr number \t Isoform ID \t TSS \t Peak start \t Peak end \t Normalized read counts
#to be counted as expressed, normalized read counts must be greater than 10 normalized reads
#to run script: python3 Pull.Norm.RNA.Counts.py <observed TSS file> <normalized read counts csv> <individual being examined: L1F, L2F, L1M, L2M, O1, O2, T1, T2> <isoform to gene id conversion file> <output file with TSSs with peaks with expression> <output file with TSSs with peaks with no expression>
#Author: Alice Naftaly, Sept 2020

import sys

#read in tss observed file
#returns dictionary with key == isoform id and value == [chr number, peak start, peak end, start_site]
def read_tss_with_peaks():
    tss_file = sys.argv[1]
    peaks_at_tss = {}
    with open(tss_file, 'r') as tss_pos:
        for line in tss_pos:
            if line.startswith("Chr"):
                continue
            else:
                new_line = line.split()
                chr_num = new_line[0]
                peak_start = int(new_line[1])
                peak_end = int(new_line[2])
                start_site = int(new_line[3])
                isoform_id = new_line[4]
                dict_value = [chr_num, peak_start, peak_end, start_site]
                if isoform_id in peaks_at_tss:
                    peaks_at_tss[isoform_id].append(dict_value)
                elif isoform_id not in peaks_at_tss:
                    peaks_at_tss.update({isoform_id:[dict_value]})
    return peaks_at_tss


#read normalized read counts file
#filters out normalized counts for a single individual
#returns dictionary with key == gene id and value == normalized read counts
def read_normalized_counts():
    counts_file = sys.argv[2]
    individual = sys.argv[3]
    all_normalized_counts = []
    gene_ids = []
    individual_counts_dict = {}
    with open(counts_file, 'r') as counts:
        for line in counts:
            new_line = line.strip("\n")
            new_line = new_line.split(",")
            if new_line[0] == "Gene_ID":
                header = new_line
            else:
                gene_id = new_line[0].strip("\"")
                gene_ids.append(gene_id)
                normalized_counts = new_line[1:len(new_line)]
                stripped_normalized_counts = [i.strip("\"") for i in normalized_counts]
                all_normalized_counts += stripped_normalized_counts
        individual_index = header.index(individual) - 1
        individual_normalized_counts = all_normalized_counts[individual_index::8]
        for ind, val in enumerate(gene_ids):
            individual_counts_dict.update({val:individual_normalized_counts[ind]})
    return individual_counts_dict

#read isoform id to gene id conversion file
#returns dictionary with key == gene id and value == isoform ids (can be mor ethan 1)
def read_conversion():
    conversion_file = sys.argv[4]
    genes_to_isoforms = {}
    with open(conversion_file, 'r') as conversion:
        for line in conversion:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                gene_id = new_line[1]
                genes_to_isoforms.update({isoform_id:gene_id})
    return genes_to_isoforms

#pull normalized reads at TSSs with peaks
#returns dictionary with key == isoform id and value == [chr number, peak start, peak end, tss, normalized read counts]
#note: isoforms from the same gene will have the same read counts (will need to use un-normalized counts or exon level analysis if need to separate these)
#I don't think it matters too much as the isoforms will likely be in the same bins and if it's the same peak it will still be counted as expressed
def pull_coverage_at_TSSs():
    genes_to_isoforms = read_conversion()
    normalized_counts = read_normalized_counts()
    TSSs = read_tss_with_peaks()
    peaks_with_expressed_genes = {}
    peaks_with_no_expression = {}
    for isoform in TSSs:
        single_gene = genes_to_isoforms[isoform]
        if single_gene in normalized_counts:
            single_gene_normalized_counts = round(float(normalized_counts[single_gene]),2)
        elif single_gene not in normalized_counts:
            single_gene_normalized_counts = 0
        single_isoform = TSSs[isoform]
        if len(single_isoform) == 1:
            single = single_isoform[0]
            if single_gene_normalized_counts >= 10:
                dict_value = [single[0], single[1], single[2], single[3], single_gene_normalized_counts]
                if isoform in peaks_with_expressed_genes:
                    peaks_with_expressed_genes[isoform].append(dict_value)
                elif isoform not in peaks_with_expressed_genes:
                    peaks_with_expressed_genes.update({isoform:[dict_value]})
            elif single_gene_normalized_counts < 10:
                dict_value = [single[0], single[1], single[2], single[3], single_gene_normalized_counts]
                if isoform in peaks_with_no_expression:
                    peaks_with_no_expression[isoform].append(dict_value)
                elif isoform not in peaks_with_no_expression:
                    peaks_with_no_expression.update({isoform:[dict_value]})
        elif len(single_isoform) > 1:
            for single in single_isoform:
                if single_gene_normalized_counts >= 10:
                    dict_value = [single[0], single[1], single[2], single[3], single_gene_normalized_counts]
                    if isoform in peaks_with_expressed_genes:
                        peaks_with_expressed_genes[isoform].append(dict_value)
                    elif isoform not in peaks_with_expressed_genes:
                        peaks_with_expressed_genes.update({isoform:[dict_value]})
                elif single_gene_normalized_counts < 10:
                    dict_value = [single[0], single[1], single[2], single[3], single_gene_normalized_counts]
                    if isoform in peaks_with_no_expression:
                        peaks_with_no_expression[isoform].append(dict_value)
                    elif isoform not in peaks_with_no_expression:
                        peaks_with_no_expression.update({isoform:[dict_value]})
    return peaks_with_expressed_genes, peaks_with_no_expression

#write output file
def write_expression():
    expressed_coverage, no_expression_coverage = pull_coverage_at_TSSs()
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = "Chr.Num\tIsoform.ID\tTSS\tPeak.Start\tPeak.End\tNormalized.Read.Counts\n"
        out.write(header)
        for isoform in expressed_coverage:
            single_isoform = expressed_coverage[isoform]
            if len(single_isoform) == 1:
                single = single_isoform[0]
                final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (single[0], isoform, single[3], single[1], single[2], single[4])
                out.write(final)
            elif len(single_isoform) > 1:
                for single in single_isoform:
                    final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (single[0], isoform, single[3], single[1], single[2], single[4])
                    out.write(final)


def write_no_expression():
    expressed_coverage, no_expression_coverage = pull_coverage_at_TSSs()
    output = sys.argv[6]
    with open(output, 'a') as out:
        header = "Chr.Num\tIsoform.ID\tTSS\tPeak.Start\tPeak.End\tNormalized.Read.Counts\n"
        out.write(header)
        for isoform in no_expression_coverage:
            single_isoform = no_expression_coverage[isoform]
            if len(single_isoform) == 1:
                single = single_isoform[0]
                final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (single[0], isoform, single[3], single[1], single[2], single[4])
                out.write(final)
            elif len(single_isoform) > 1:
                for single in single_isoform:
                    final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (single[0], isoform, single[3], single[1], single[2], single[4])
                    out.write(final)


#call all functions
def call():
    expressed_TSSs = write_expression()
    no_expression = write_no_expression()

call()
