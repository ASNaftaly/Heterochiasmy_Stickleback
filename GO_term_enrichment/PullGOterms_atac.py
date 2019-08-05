#Pulling GO terms from csv file
#pulls GO terms for all genes; if no GO terms available, will have n/a
#to run script: python3 PullGOterms_atac.py <csv file> <output file in tab delimited for so it can be read into excel file>
#Author: Alice Naftaly, July 2019

import sys

#pulls GO terms for each gene
#this won't pull genes that have no GO terms
#returns dictionary with key = Gene ID and value = GO terms
def pull_GO_ids():
    csv_file = sys.argv[1]
    gene_dict = {}
    gene_id = ""
    GO_term = ""
    x = 0
    with open(csv_file, 'r') as csv:
        for line in csv:
            new_line = line.split(",")
            for value in new_line:
                if value.startswith("ENSGACG"):
                    gene_id = value
                elif value.startswith("GO:"):
                    GO_term = value
                    x += 1
            if gene_id in gene_dict:
                gene_dict[gene_id].append(GO_term)
            elif gene_id not in gene_dict:
                gene_dict.update({gene_id:[GO_term]})
    return gene_dict

#pulling genes that have no GO terms
#returns dictionary with key = Gene ID and value = "n/a"
#in case these are needed for something else
def no_GO_terms():
    csv_file = sys.argv[1]
    gene_dict = {}
    gene_id = ""
    with open(csv_file, 'r') as csv:
        for line in csv:
            x = 0
            new_line = line.split(",")
            if new_line[0].startswith("Gene"):
                continue
            else:
                for value in new_line:
                    if value.startswith("ENSGACG"):
                        gene_id = value
                    elif value.startswith("GO:"):
                        x = 1
                if x == 0:
                    if gene_id in gene_dict:
                        gene_dict[gene_id].append("n/a")
                    elif gene_id not in gene_dict:
                        gene_dict.update({gene_id:["n/a"]})
            x = 0
    return gene_dict

#writing final dictionary to tab delimited file
#format = GENE ID \t GO terms in list (separated by commas)
def write():
    final_dict = pull_GO_ids()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for key in final_dict:
            single_value = final_dict[key]
            final_value = ",".join(single_value)
            final = "%s\t%s\n" % (str(key), str(final_value))
            out.write(final)

write()
