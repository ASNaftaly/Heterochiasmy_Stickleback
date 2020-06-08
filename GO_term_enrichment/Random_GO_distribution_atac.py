#Compare GO terms for genes at peaks
#1. randomly pull specific number of genes from all genes list while excluding the genes that are in the actual data set that the random distribution is being created for
#2. pull those GO terms and return output with summary of GO terms found
#   this will have GO term \t number of occurrences in permutations
#3. repeat 10,000 times
#   this file will be added to for all 10,000 permutations
#   then GO terms can be expanded to 10,000 occurrences where if the GO term was not counted in a permutation that will be added as 0 because each permutation will only have the GO term recorded once in summary file with the total number of times it was seen in that permutation
#4. Compare to known GO terms (somehow; will explain in another script)
#to run script: python3 Random_GO_distribution_atac.py <all genes GO term file with Gene IDs> <actual test file in format of GO terms *has Gene ID> <random distribution output>
#Author: Alice Naftaly, July 2019

import sys
import random

#read all genes file
#format: Gene ID \t GO terms in list (separated by commas)
def pull_genes():
    all_genes_file = sys.argv[1]
    gene_dict = {}
    with open(all_genes_file, 'r') as all_genes:
        for line in all_genes:
            new_line = line.split("\t")
            gene_ID = new_line[0]
            GO_term = new_line[1].strip("\n")
            split_go_terms = GO_term.split(",")
            if gene_ID in gene_dict:
                gene_dict[gene_ID].append(split_go_terms)
            elif gene_ID not in gene_dict:
                gene_dict.update({gene_ID:split_go_terms})
    return gene_dict

#pull genes from observed set
#returns list of observed gene Ids
#number of observed genes should equal the number of genes near peaks that had GO terms
def pull_observed_genes():
    observed_file = sys.argv[2]
    gene_ids = []
    with open(observed_file, 'r') as observed:
        for line in observed:
            new_line = line.split()
            gene_ids.append(new_line[0])
    return gene_ids

#returns number of genes in observed set and the number to randomly assign for each permutation
#number of observed genes should equal the number of genes near peaks that had GO terms
#don't need to include the ones that have no GO terms (may revist this by including these and seeing how this affects the results)
def num_permutations():
    observed_genes = pull_observed_genes()
    return len(observed_genes)

#filter all genes list to remove genes in observed gene set
#returns all genes that are not part of observed set
def filter_all_genes():
    observed_genes = pull_observed_genes()
    all_genes = pull_genes()
    for value in observed_genes:
        if value in all_genes.keys():
            del all_genes[value]
    return all_genes

#randomly select a subset of genes from gene_dict
#returns list of GO terms found in random subset of genes
def random_genes():
    genes = pull_genes()
    gene_names = list(genes.keys())
    total_genes = num_permutations()
    random_go_terms = []
    random_go_terms_dict = {}
    x = 0
    while x < int(total_genes):
        random_gene = random.choice(gene_names)
        random_go_terms += genes[random_gene]
        for value in random_go_terms:
            key = value
            dict_value = "1"
            if key in random_go_terms_dict:
                random_go_terms_dict[key].append(dict_value)
            elif key not in random_go_terms_dict:
                random_go_terms_dict.update({key:[dict_value]})
        x += 1
    return random_go_terms_dict

#write all GO terms to one file
def write():
    random_go_terms_dict = random_genes()
    output_file = sys.argv[3]
    with open(output_file, 'a') as output:
        for key in random_go_terms_dict:
            go_term = key
            num_occurrences = len(random_go_terms_dict[key])
            #final output = GO term \t total number of occurrences for single permutations
            final = "%s\t%s\n" % (str(go_term),str(num_occurrences))
            output.write(final)

write()
