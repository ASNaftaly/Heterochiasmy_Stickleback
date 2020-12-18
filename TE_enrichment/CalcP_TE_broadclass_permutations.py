#calculate p values for TE enrichment (broad classes; DNA, LINE, SINE, LTR, and RC)
#will need to read in permutations file with 1 total per permutation (combined all autosomes in ExpectedDist_atac_TEbroadcounts.py)
#want to know how many permutations are smaller than the observed TE enrichment
#to run script: CalcP_TE_broadclass_permutations.py <permutations file from ExpectedDist_atac_TEbroadcounts.py> <TE type: DNA, LINE, SINE, LTR, RC> <observed TE count> <output>
#Author: Alice Naftaly, Oct 2020

import sys

#pull random summary
#each line = 1 permutation with total TEs covered
#returns list of permutation totals
def pull_random_summary():
	permutations_file = sys.argv[1]
	permutations_list = []
	with open(permutations_file, 'r') as permutations:
		for line in permutations:
			new_line = int(line.strip())
			permutations_list.append(new_line)
	return permutations_list

#comparing observed values with random permutations
#returns pvalue
def compare():
	permutations = pull_random_summary()
	observed = int(sys.argv[3])
	perms_below_observed = 0
	for perm in permutations:
		if perm < observed:
			perms_below_observed += 1
		pvalue = round(1 - (perms_below_observed/10000), 6)
	return pvalue

#write p values to file
def write():
	pvalue = compare()
	te_type = sys.argv[2]
	output = sys.argv[4]
	with open(output, 'a') as out:
		final = "%s\t%s\n" % (str(te_type), str(pvalue))
		out.write(final)


write()
