#randomize hotspots to figure out TSS enrichment using new TSS positions
#will need to read in chromosome sizes
#to run script: python3 Randomize.Hotspots.py <Chromosome sizes file> <file with num of hotspots to randomize per chr> <output file: chr num random hot start random hot end> <output file>
#Author: Alice Naftaly, Oct 2020

import sys
import random

#read in chromosome sizes
#returns dictionary with key == chr number and value == chr size
def read_chr_sizes():
    size_file = sys.argv[1]
    size_dict = {}
    with open(size_file, 'r') as sizes:
        for line in sizes:
            new_line = line.split()
            chr_num = new_line[0]
            chr_size = int(new_line[1])
            size_dict.update({chr_num:chr_size})
    return size_dict


#read in number of hotspots to randomize
#returns dictionary with key = chr number and value == number of hotspots to randomize
def pull_randomization_number():
    hot_file = sys.argv[2]
    numbers_dict = {}
    with open(hot_file, 'r') as hot:
        for line in hot:
            new_line = line.split()
            chr_num = new_line[0]
            num_to_randomize = int(new_line[1])
            numbers_dict.update({chr_num:num_to_randomize})
    return numbers_dict


#simulates random spots based on size of chromosome and the number of spots identified from each chromosome
def randomize_spots():
    chr_sizes = read_chr_sizes()
    random_nums = pull_randomization_number()
    random_list = []
    random_spots = {}
    for key in random_nums:
        single_chr_size = chr_sizes[key]
        single_num_to_randomize = random_nums[key]
        x = 0
        while x < single_num_to_randomize:
            random_spot = random.randrange(0,single_chr_size)
            if random_spot in random_list:
                continue
            else:
                random_list.append(int(random_spot))
                random_spot_start = random_spot
                random_spot_end = int(random_spot) + 2000
                single_spot_value = [random_spot_start, random_spot_end]
            if key in random_spots:
                random_spots[key].append(single_spot_value)
            elif key not in random_spots:
                random_spots.update({key:[single_spot_value]})
            x += 1
        random_list = []
    return random_spots


#write random peaks to file
#output = chr num \t random peak start \t random peak end \n
def write():
    random_spots = randomize_spots()
    output_file = sys.argv[3]
    with open(output_file, 'w') as out:
        for key in random_spots:
            single_value = random_spots[key]
            for value in single_value:
                start = value[0]
                end = value[1]
                final = "%s\t%s\t%s\n" % (str(key), str(start), str(end))
                out.write(final)


write()
