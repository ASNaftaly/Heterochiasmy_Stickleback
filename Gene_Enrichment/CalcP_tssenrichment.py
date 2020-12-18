#calculating p value for enrichment of random peaks at TSSs
#to run script: python3 CalcP_tssenrichment.py < output from Total_permutations_autosomes_TSS_enrichment.py> <observed peaks in TSSs; integer> <sample: L1F, L2F, L1M, L2M, O1, O2, T1, T2> >> <output to standard out>
#Author: Alice Naftaly, Oct 2020

import sys

#pull random summary
#each line = 1 total random peaks per permutation
#returns list with total random peaks
def pull_random_summary():
	permutations_file = sys.argv[1]
	permutations_list = []
	with open(permutations_file, 'r') as permutations:
		for line in permutations:
			new_line = int(line.strip("\n"))
			permutations_list.append(new_line)
	return permutations_list


#comparing observed values with random permutations
def compare():
	permutations = pull_random_summary()
	observed = int(sys.argv[2])
	sample = sys.argv[3]
	perms_below_observed = 0
	for p in permutations:
		if p <= observed:
			perms_below_observed += 1
	pvalue = round(1 - (perms_below_observed/10000), 6)
	final = "%s p-value: %s\n" % (str(sample), str(pvalue))
	print(final)


compare()
