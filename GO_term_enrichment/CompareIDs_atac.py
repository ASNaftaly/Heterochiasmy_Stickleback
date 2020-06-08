#Comparing Gene IDs pulled from TSSs near peaks
#to run script: python3 CompareIDs.py <id file 1> <id file 2> <output file 1, shared genes> <output file 2, unique file 1 genes> <output file 3, unique file 2 genes>
#Author: Alice Naftaly, July 2019, edited from CompareIDs.py

import sys

#pulling IDS from file 1
def pull_ids_1():
    id_file = sys.argv[1]
    id_dict = {}
    with open(id_file, 'r') as ids:
        for line in ids:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                final_id = new_line[2]
                if chr_num in id_dict:
                    id_dict[chr_num].append(final_id)
                elif chr_num not in id_dict:
                    id_dict.update({chr_num:[final_id]})
    return id_dict

#pulling IDS from file 2
def pull_ids_2():
    id_file = sys.argv[2]
    id_dict = {}
    with open(id_file, 'r') as ids:
        for line in ids:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                final_id = new_line[2]
                if chr_num in id_dict:
                    id_dict[chr_num].append(final_id)
                elif chr_num not in id_dict:
                    id_dict.update({chr_num:[final_id]})
    return id_dict

#comparing IDs from files
def compare():
    ids_1 = pull_ids_1()
    ids_2 = pull_ids_2()
    output_shared = sys.argv[3]
    output_unique_1 = sys.argv[4]
    output_unique_2 = sys.argv[5]
    with open(output_shared, 'a') as out, open(output_unique_1, 'a') as out2, open(output_unique_2, 'a') as out3:
        for key in ids_1:
            if key in ids_2:
                id_list_1 = ids_1[key]
                id_list_2 = ids_2[key]
                for value in id_list_1:
                    if value in id_list_2:
                        final_shared = "%s\n" % (str(value))
                        out.write(final_shared)
                    elif value not in id_list_2:
                        final_ids_1 = "%s\n" % (str(value))
                        out2.write(final_ids_1)
                for val in id_list_2:
                    if val not in id_list_1:
                        final_ids_2 = "%s\n" % (str(val))
                        out3.write(final_ids_2)

compare()
