#pull TSSs that have no peak (within 500bp) but are expressed
#need to read in observed peak TSS overlap file, all TSSs
#output file format:
#chr number \t Isoform ID \t TSS \t Peak start \t Peak end \t Normalized read counts
#to run script: python3 Pull.TSSs.with.No.Peaks.py <observed TSS file> <normalized read counts csv> <individual being examined: L1F, L2F, L1M, L2M, O1, O2, T1, T2> <isoform to gene id conversion file> <all isoforms bed> <output file>

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
                peaks_at_tss.update({isoform_id:dict_value})
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
                if gene_id in genes_to_isoforms:
                    genes_to_isoforms[gene_id].append(isoform_id)
                elif gene_id not in genes_to_isoforms:
                    genes_to_isoforms.update({gene_id:[isoform_id]})
    return genes_to_isoforms

#all isoforms bed file
#Need isoform id, chr number, and TSS
def read_all_isoforms_bed():
    bed_file = sys.argv[5]
    isoforms_dict = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            new_line = line.split()
            chr_num = new_line[0]
            isoform_id = new_line[3]
            strand = new_line[5]
            if strand == "+":
                tss = new_line[1]
            elif strand == "-":
                tss = new_line[2]
            dict_value = [chr_num, tss]
            isoforms_dict.update({isoform_id:dict_value})
    return isoforms_dict

#pull TSSs with no peaks
#returns dictionary with key == isoform id and value = []
def pull_TSSs_with_no_peaks():
    genes_to_isoforms = read_conversion()
    normalized_counts = read_normalized_counts()
    TSSs = read_tss_with_peaks()
    bed_dict = read_all_isoforms_bed()
    no_peaks_with_expression = {}
    for gene in genes_to_isoforms:
        single_gene_isoforms = genes_to_isoforms[gene]
        if gene in normalized_counts:
            single_gene_normalized_counts = round(float(normalized_counts[gene]),2)
            for isoform in single_gene_isoforms:
                if isoform not in TSSs and single_gene_normalized_counts > 10:
                    if isoform in bed_dict:
                        single_bed = bed_dict[isoform]
                        dict_value = [single_bed[0], single_bed[1], single_gene_normalized_counts]
                        no_peaks_with_expression.update({isoform:dict_value})
    return no_peaks_with_expression


#write output file
def write():
    no_peaks_with_expression = pull_TSSs_with_no_peaks()
    output = sys.argv[6]
    with open(output, 'a') as out:
        header = "Chr.Num\tIsoform.ID\tTSS\tNormalized.Read.Counts\n"
        out.write(header)
        for isoform in no_peaks_with_expression:
            single_isoform = no_peaks_with_expression[isoform]
            final = "%s\t%s\t%s\t%s\n" % (single_isoform[0], isoform, single_isoform[1], single_isoform[2])
            out.write(final)


write()
