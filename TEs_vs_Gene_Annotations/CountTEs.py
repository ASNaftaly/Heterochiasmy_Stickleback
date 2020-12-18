#counting TEs
#to run script: python3 CountTEs.py <repeat masker output>
#Author: Alice Naftaly, Sept 2020

import sys

#pull class identifiers
def pull_TEclasses():
    te_file = sys.argv[1]
    te_counts = {}
    with open(te_file, 'r') as tes:
        for line in tes:
            new_line = line.split()
            if len(new_line) == 0 or new_line[0] == "SW" or new_line[0] == "score" or new_line[(len(new_line)-1)] == "*":
                continue
            else:
                te_class = new_line[10]
                if te_class in te_counts:
                    te_counts[te_class].append("1")
                elif te_class not in te_counts:
                    te_counts.update({te_class:["1"]})
    count = 0
    for te in te_counts:
        print(te)
        print(len(te_counts[te]))
        count += len(te_counts[te])
    print(count)

pull_TEclasses()
