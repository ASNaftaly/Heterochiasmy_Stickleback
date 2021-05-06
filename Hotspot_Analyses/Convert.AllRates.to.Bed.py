#changing allrates file from hotspot paper to bed file format to convert positions
#all rates file format: Position \t Recomb.rate
#need the following format: chr \t start.pos \t end.pos \t recombination rates
#where start pos is the first position that has a specific recombination rate and the end pos is the last position (in sequence!) that has the same rate
#to run script: python3 Convert.AllRates.to.Bed.py <all rates file from a single chromosome> <chromosome number as chr roman numeral> <output file>

import sys
from collections import OrderedDict


#read all rates file
#returns dictionary with key == recombination rate and value == list of positions with that recombination rate
def read_allrates():
    rates_file = sys.argv[1]
    rates_dict = {}
    with open(rates_file, 'r') as rates:
        for line in rates:
            if line.startswith("Position"):
                continue
            else:
                new_line = line.split()
                pos = int(new_line[0])
                recomb_rate = new_line[1]
                if recomb_rate in rates_dict:
                    rates_dict[recomb_rate].append(pos)
                elif recomb_rate not in rates_dict:
                    rates_dict.update({recomb_rate:[pos]})
    return rates_dict

#pull first and last position with the same recombination rate
#if same recombination rate shows up in multiple places in the genome will need to have multiple start and end positions
def condense_positions():
    rates = read_allrates()
    condensed_rates_dict = {}
    for rate in rates:
        single_rate = rates[rate]
        x = 0
        condensed_dict = {}
        count = 1
        while x < len(single_rate)-1:
            y = single_rate[x + 1]
            if y == single_rate[x] + 1:
                if count in condensed_dict:
                    condensed_dict[count].append(single_rate[x])
                elif count not in condensed_dict:
                    condensed_dict.update({count:[single_rate[x]]})
                x += 1
            elif y != single_rate[x] + 1:
                count += 1
                condensed_dict.update({count:[single_rate[x+1]]})
                x += 1
        for key in condensed_dict:
            single_key = condensed_dict[key]
            dict_value = [min(single_key), max(single_key)+1]
            if rate in condensed_rates_dict:
                condensed_rates_dict[rate].append(dict_value)
            elif rate not in condensed_rates_dict:
                condensed_rates_dict.update({rate:[dict_value]})
    return condensed_rates_dict


#ordering rates by positions
#returns dictionary with key == position (start position for region with the same recombination rate) and value == [start position, end position, rate]
def order_positions():
    condensed_dict = condense_positions()
    positions_dict = {}
    for key in condensed_dict:
        single_region = condensed_dict[key]
        if len(single_region) == 1:
            start = single_region[0][0]
            end = single_region[0][1]
            dict_value = [start, end, key]
            positions_dict.update({start:dict_value})
        elif len(single_region) > 1:
            for single in single_region:
                start = single[0]
                end = single[1]
                dict_value = [start, end, key]
                positions_dict.update({start:dict_value})
    ordered_dict = OrderedDict(sorted(positions_dict.items(),key=lambda t:t[0]))
    return ordered_dict

#write output file
def write():
    final_data = order_positions()
    chr_num = sys.argv[2]
    output_file = sys.argv[3]
    with open(output_file, 'a') as out:
        for key in final_data:
            single_region_data = final_data[key]
            final = "%s\t%s\t%s\t%s\n" % (chr_num, str(single_region_data[0]), str(single_region_data[1]), single_region_data[2])
            out.write(final)

write()
