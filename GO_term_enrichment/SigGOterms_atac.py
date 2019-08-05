#Determining if GO terms are enriched in peaks near TSSs
#1. Read random permutations file which is the format:
#   GO term \t 10,000 occurrences ranging from 0 and up
#2. read in observed go term occurrences for specfic observed file and save as dictionary with key = GO term and value = number of observed occurrences
#5. Go key by key and ask how many times the observed value is greater then the random value; get total and use 1- (total/10000) to get pvalue
#6. Write final file: GO term \t pvalue
#might need to adjust sig p value for each set depending on number of tests completed
#can use bonferroni correction?
#to run script: python3 SigGoTerms_atac.py <full random dist GO terms> <observed GO terms file> <final output file>
#Author: Alice Naftaly, July 2019

import sys

#creates dictionary from GO terms from random permutations
#key = GO term
#value = 1 = how many times GO term was present (pull length of key to get total counts of each GO term)
def pull_terms():
    go_terms = "L1.unique.random.dist.GO.terms.10k.permutations.summary.txt"#sys.argv[1]
    go_dict = {}
    with open(go_terms, "r") as go:
        for line in go:
            new_line = line.split("\t")
            key = new_line[0]
            occurrences = new_line[1:len(new_line)-1]
            for value in occurrences:
                if key in go_dict:
                    go_dict[key].append(int(value))
                elif key not in go_dict:
                    go_dict.update({key:[int(value)]})
    return go_dict

#creates dictionary from observed GO term occurrences (counts file)
#key = GO term
#value = 1 = how many times GO term was present (pull length of key to get total counts of each GO term)
def pull_observed_terms():
    go_terms = "L1.sigpeaks.unique.Ensembl97.GO.term.counts.txt"#sys.argv[2]
    obs_go_dict = {}
    with open(go_terms, 'r') as obs_go:
        for line in obs_go:
            if line.startswith("GO:"):
                new_line = line.split("\t")
                key = new_line[0]
                value = new_line[1].strip("\n")
                obs_go_dict.update({key:value})
    return obs_go_dict

#comparing files
#calculates p value for each GO term
#just need to correct pvalue for multiple tests
def compare():
    observed = pull_observed_terms()
    random_dist = pull_terms()
    single_go_term_count = 0
    output = "Test.txt"#sys.argv[3]
    with open(output, 'a') as out:
        header = "GO.Term\tP.Value\tPresentinRandDist\n"
        out.write(header)
        for key in observed:
            if key in random_dist:
                #print("GO term")
                #print(key)
                #total occurrences in random distribution
                total_occur_random = random_dist[key]
                for value in total_occur_random:
                    if int(value) <= int(observed[key]):
                        single_go_term_count += 1
                if single_go_term_count == 0:
                    print(key)
                #print("Number of permutations counts that are smaller than observed")
                #print(single_go_term_count)
                pvalue = round(1 - (single_go_term_count/10000),5)
                final = "%s\t%s\tyes\n" % (str(key),str(pvalue))
                #out.write(final)
                single_go_term_count = 0
            else:
                final = "%s\t%s\tno\n" % (str(key),str(0))
                out.write(final)

compare()
