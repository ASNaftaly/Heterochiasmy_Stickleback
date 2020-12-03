#pulling top 20 upregulated testis genes compared to ovary samples
#got normalized counts from DESeq2 in R
#got adjusted p values and log2fold changes from DESeq2 in R
#testis upregulated genes will have a negative log2fold change as these changes are in terms of ovary expression
#normalized counts and adjusted pvalues are in csv files
#to run script: python3 Pull_top_DE_genes.py <normalized read counts csv> <padj values in csv> <output file>
#Author: Alice Naftaly, Dec 2020

import sys

#read in normalized counts
#returns dictionary with key == gene and value == [average ovary normalized read count, average testis normalized read count]
def read_normalized_counts():
    normalized_file = sys.argv[1]
    normalized_counts_dict = {}
    with open(normalized_file, 'r') as normalized_counts:
        for line in normalized_counts:
            new_line = line.split(",")
            if new_line[0].startswith("Genes"):
                continue
            else:
                gene = new_line[0].strip("\"")
                #values in list are in the following order:
                #L1F, L2F, L1M, L2M, O1, O2, T1, T2
                #want average of Ovary and testis
                ave_ovary_norm_read_counts = (float(new_line[5]) + float(new_line[6]))/2
                ave_testis_norm_read_counts = (float(new_line[7]) + float(new_line[8].strip("\n")))/ 2
                ave_norm_values = [ave_ovary_norm_read_counts, ave_testis_norm_read_counts]
                normalized_counts_dict.update({gene:ave_norm_values})
    return normalized_counts_dict

#read in adjusted p values and log2fold values
#returns dictionary where key == gene and value == [log2_fold, adjusted p value]
def read_padj_values():
    padj_file = sys.argv[2]
    padj_dict = {}
    with open(padj_file, 'r') as padj_vals:
        for line in padj_vals:
            new_line = line.split(",")
            padj = ""
            if new_line[0] == "genes":
                continue
            else:
                gene = new_line[0].strip("\"")
                log2_fold = float(new_line[1])
                if new_line[2].strip("\n") == "NA":
                    continue
                else:
                    padj = float(new_line[2].strip("\n"))
                dict_value = [log2_fold, padj]
                padj_dict.update({gene:dict_value})
    return padj_dict


#sorting padj dict to include the top genes with padj < 10e-6 and have negative log2_fold (meaning greater expression in the testis)
#returns dictionary with key == gene and value == [log2fold value, padj]
def sort_padj():
    padj_dict = read_padj_values()
    significant_padj = {}
    for gene in padj_dict:
        single_gene = padj_dict[gene]
        #selects for genes with adjusted pvalue below threshold
        if single_gene[1] < 10e-6:
            #selects genes that have negative log2fold values
            if single_gene[0] < 0:
                significant_padj.update({gene:single_gene})
    return significant_padj

#write output
def write():
    sorted_padj = sort_padj()
    normalized_reads = read_normalized_counts()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Gene\tPadj\tAve.Ovary.Norm.Read.Count\tAve.Testis.Norm.Read.Count\n"
        out.write(header)
        for gene in sorted_padj:
            single_padj_gene = sorted_padj[gene]
            adjust_pvalue = single_padj_gene[1]
            normalized_read_counts = normalized_reads[gene]
            ovary_norm = normalized_read_counts[0]
            testis_norm = normalized_read_counts[1]
            final = "%s\t%s\t%s\t%s\n" % (str(gene), str(adjust_pvalue), str(ovary_norm), str(testis_norm))
            out.write(final)

write()
