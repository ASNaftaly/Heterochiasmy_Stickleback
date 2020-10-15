#pull hotspot start position from final hot file
#format: separated per chromosome and where each line == Hotspot Window: start position
#need to pull start position and add 2kb for end position
#should add chromosome to line as well
#final output =
#chr number \t hot start \t hot end
#to run script: python3 PullHotPos.py <hot file> <chr number in roman numerals> <output file>
#Author: Alice Naftaly, Oct 2020

import sys

#read in hot file
def read_hot():
    hot_file = sys.argv[1]
    hotspot_list = []
    with open(hot_file, 'r') as hot_spots:
        for line in hot_spots:
            if line.startswith("Hotspot"):
                new_line = line.split()
                hot_start = int(new_line[2])
                hot_end = hot_start + 2000
                hot = [hot_start, hot_end]
                hotspot_list.append(hot)
    return hotspot_list


#write hotspots
def write():
    hotspots = read_hot()
    chr_num = sys.argv[2]
    output = sys.argv[3]
    with open(output, 'a') as out:
        for hot in hotspots:
            hot_start = str(hot[0])
            hot_end = str(hot[1])
            final = "%s\t%s\t%s\n" % (str(chr_num), hot_start, hot_end)
            out.write(final)

write()
