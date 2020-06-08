#Counting the number of occurrences for each GO term in 10,000 permutations
#Each time a GO term is present in the file, this shows how many times out of 10,000 that this GO term was found
#values after each GO term represent how many times that GO term was found in a given permutation
#writes values to a output file as follows:
#   GO term \t #occurrences permutation 1 \t #occurrences permutation 2 \t ...
#to run script: python3 Summarize_GO_permutations.py <random dist of GO occurrences> <output summary file>
#Author: Alice Naftaly, July 2019

import sys

#pulls go terms and returns dictionary
#dictionary has key = go term and value = number of times that GO term is seen in a given permutation
def pull_GO_terms():
    permutations_file = sys.argv[1]
    go_dict = {}
    with open(permutations_file, 'r') as permutations:
        for line in permutations:
            new_line = line.split()
            key = new_line[0]
            value = new_line[1]
            if key in go_dict:
                go_dict[key].append(value)
            elif key not in go_dict:
                go_dict.update({key:[value]})
    return go_dict

#expand number of values for in key to 10,000
def expand():
    go_dict = pull_GO_terms()
    for key in go_dict:
        length_start = len(go_dict[key])
        while length_start <= 10000:
            go_dict[key].append("0")
            length_start += 1
    return go_dict

#writing output to file that can be easily read for determining significance
def write():
    go_dict = expand()
    output_file = sys.argv[2]
    with open(output_file, 'a') as output:
        for key in go_dict:
            occurrences = go_dict[key]
            string_occurrences = "\t".join(occurrences)
            final = str(key) + "\t" + string_occurrences + "\n"
            output.write(final)

write()
