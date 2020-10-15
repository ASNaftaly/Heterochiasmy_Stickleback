#pulling sequence at hotspots to realign hotspots to V5 assembly for a few analyses
#position format:
#chr num \t hot start \t hot end
#need to read in one chromosome at a time, stickleback fasta, and hotspot positions
#will only examine autosomes
#to run script: python3 Pull.Sequences.at.Hotspots.py <hotspots position file> <chr number> <stickleback fasta>  <output fasta file>
#Author: Alice Naftaly, October 2020


import sys

#read in hot positions file:
def read_hot():
    hot_file = sys.argv[1]
    chr = sys.argv[2]
    hot_list = []
    with open(hot_file, 'r') as hotspots:
        for line in hotspots:
            if line.startswith("chr"):
                new_line = line.split()
                chr_num = new_line[0]
                if chr == chr_num:
                    hot_start = int(new_line[1])
                    hot_end = int(new_line[2])
                    hot_spot = [hot_start, hot_end]
                    hot_list.append(hot_spot)
    return hot_list


#read stickleback fasta
def read_fasta():
    fasta_file = sys.argv[3]
    chr = sys.argv[2]
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.strip('>')
                final_fasta_id = new_line.strip("\n")
            else:
                new_line = line.strip("\n")
                if final_fasta_id in fasta_dict:
                    fasta_dict[final_fasta_id].append(new_line)
                elif final_fasta_id not in fasta_dict:
                    fasta_dict.update({final_fasta_id:[new_line]})
    fasta_seq = []
    for l in fasta_dict[chr]:
        fasta_seq += l
    return fasta_seq

#pull sequences at peaks:
#returns dictionary with key == peak.startpos.endpos and value == sequence at peak
def pull_seqs():
    hotspots = read_hot()
    chr_num = sys.argv[2]
    single_chr_fasta = read_fasta()
    hot_sequences = {}
    for index, hot in enumerate(hotspots):
        single_hot_header = "Hotspot.%s.%s" % (str(chr_num), str(index))
        hot_start = hot[0]
        hot_end = hot[1]
        single_hot_sequence = single_chr_fasta[hot_start:hot_end]
        hot_sequences.update({single_hot_header:single_hot_sequence})
    return hot_sequences


#write output fasta
def write():
    hot_seqs = pull_seqs()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for hot in hot_seqs:
            header = ">" + hot + "\n"
            sequence = hot_seqs[hot]
            joined_seq = "".join(sequence)
            out.write(header)
            out.write(joined_seq + "\n")

write()
