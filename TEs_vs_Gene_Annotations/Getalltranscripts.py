#script to pull genes not found in isoseq dataset (therefore not in liver, brain, pronephros, and gonads) from ensembl build 97 genes for TE/TSS overlap
#will need to read in ensembl gtf lifted to V5, isoseq classification, all isoforms bed file
#output should be in bed format
#chr_number start_pos end_pos isoform_id or transcript_id "." strand
#to run script: python3 Getalltranscripts.py <isoseq bed file with V5 coordinates> <isoseq classification file>

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


#read in classification file
#need isoform id and gene id only
#returns list with key == gene id
def read_class():
    class_file = sys.argv[2]
    isoseq_genes = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                gene_id = new_line[6]
                isoseq_genes.append(gene_id)
    return isoseq_genes


#pull genes from ensembl GTF
#need gene id, transcript id, start position, end position, and Strand
#returns dictionary with key == gene id and value == for each transcript [chr_num, transcript id, start pos, end pos, strand]
def read_ensembl_gtf():
    gtf_file = sys.argv[3]
    gtf_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            feature = new_line[2]
            if feature == "transcript":
                chr_num = new_line[0]
                start_pos = new_line[3]
                end_pos = new_line[4]
                strand = new_line[6]
                transcript_info = new_line[8].split(";")
                full_gene_id = transcript_info[0].split(" ")
                gene_id = full_gene_id[1].strip("\"")
                full_transcript_id = transcript_info[2].split(" ")
                transcript_id = full_transcript_id[2].strip("\"")
                dict_value = [chr_num, transcript_id, start_pos, end_pos, strand]
                if gene_id in gtf_dict:
                    gtf_dict[gene_id].append(dict_value)
                elif gene_id not in gtf_dict:
                    gtf_dict.update({gene_id:[dict_value]})
    return gtf_dict

#sort ensembl gtf based on genes in isoseq_genes
#returns dictionary with genes that are not found in the isoseq data set
def sort_ensembl():
    ensembl_genes = read_ensembl_gtf()
    isoseq_genes = read_class()
    final_ensembl_genes = {}
    for gene in ensembl_genes:
        if gene not in isoseq_genes:
            final_ensembl_genes.update({gene:ensembl_genes[gene]})
    return final_ensembl_genes



#write bed output
def write():
    ensembl_genes = sort_ensembl()
    isoseq_genes = read_isoform_bed()
    chrs = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrXVII", "chrXVIII", "chrXIX", "chrXX", "chrXXI", "chrY"]
    output = sys.argv[4]
    with open(output, 'a') as out:
        for chr in chrs:
            single_isoseq = isoseq_genes[chr]
            for iso_value in single_isoseq:
                final_iso = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chr), iso_value[1], iso_value[2], iso_value[0], ".", iso_value[3])
                out.write(final_iso)
        for gene in ensembl_genes:
            single_gene = ensembl_genes[gene]
            for ens_value in single_gene:
                final_ens = "%s\t%s\t%s\t%s\t%s\t%s\n" % (ens_value[0], ens_value[2], ens_value[3], ens_value[1], ".", ens_value[4])
                out.write(final_ens)
write()
