#identifying genes with normalized read count above 10 in testis samples only
#read in csv file with all genes and all samples
#to run script: python3 Identify.testis.expressed.genes.py <csv file with all normalized RNA counts>

import sys

#read in csv file
#returns dictionary with key == gene id and value == normalized read counts for all samples
def read_csv():
    csv_file = sys.argv[1]
    expression_dict = {}
    with open(csv_file, 'r') as csv:
        for line in csv:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split(",")
                gene_id = new_line[0].strip("\"")
                l1f_nrc = float(new_line[1].strip("\""))
                l2f_nrc = float(new_line[2].strip("\""))
                l1m_nrc = float(new_line[3].strip("\""))
                l2m_nrc = float(new_line[4].strip("\""))
                o1_nrc = float(new_line[5].strip("\""))
                o2_nrc = float(new_line[6].strip("\""))
                t1_nrc = float(new_line[7].strip("\""))
                t2 = new_line[8].strip("\n")
                t2_nrc = float(t2.strip("\""))
                dict_value = [l1f_nrc, l2f_nrc, l1m_nrc, l2m_nrc, o1_nrc, o2_nrc, t1_nrc, t2_nrc]
                expression_dict.update({gene_id:dict_value})
    return expression_dict


#compare normalized read counts across tissues
def find_testis_specific_genes():
    expression_dict = read_csv()
    testis_specific_genes = []
    testis_upregulated_genes = []
    for gene in expression_dict:
        single_gene = expression_dict[gene]
        l1f_nrc = single_gene[0]
        l2f_nrc = single_gene[1]
        l1m_nrc = single_gene[2]
        l2m_nrc = single_gene[3]
        o1_nrc = single_gene[4]
        o2_nrc = single_gene[5]
        t1_nrc = single_gene[6]
        t2_nrc = single_gene[7]
        if l1f_nrc < 10 and l2f_nrc < 10 and l1m_nrc < 10 and l2m_nrc < 10 and o1_nrc < 10 and o2_nrc < 10 and t1_nrc >= 10 and t2_nrc >= 10:
            testis_specific_genes.append(gene)
        elif t1_nrc > 2*l1f_nrc and t1_nrc > 2*l2f_nrc and t1_nrc > 2* l1m_nrc and t1_nrc > 2*l2m_nrc and t1_nrc > o1_nrc and t1_nrc > 2*o2_nrc and t2_nrc > 2*l1f_nrc and t2_nrc > 2*l2f_nrc and t2_nrc > 2* l1m_nrc and t2_nrc > 2*l2m_nrc and t2_nrc > o1_nrc and t2_nrc > 2*o2_nrc:
            testis_upregulated_genes.append(gene)
    return testis_specific_genes, testis_upregulated_genes

#write testis specific genes
def write_testis_specific():
    t_genes, t_upregulated_genes =find_testis_specific_genes()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for gene in t_genes:
            final = "%s\n" % str(gene)
            out.write(final)


#write testis upregulated genes
def write_testis_upregulated():
    t_genes, t_upregulated_genes =find_testis_specific_genes()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for gene in t_upregulated_genes:
            final = "%s\n" % str(gene)
            out.write(final)

#call functions
def call():
    testis_specific = write_testis_specific()
    testis_upregulated = write_testis_upregulated()

call()
