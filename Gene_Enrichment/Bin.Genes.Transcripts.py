#breaking down genes and transcripts into bins (matching bins created for accessible regions)
#each bin represents 1 twentieth of a chromosome where the centromere position is taken into account; 20 bins per chromosome
#will create 20 bins but simplify to 10 bins (representing a folded chromosome at the centromere)
#output = single file with: Chr.Num \t Bin Number \t Bin Start \t Bin End \t Gene or Transcript Name \n
#to run script: python3 Bin.Genes.Transcripts.py <chromosome sizes file: Chr \t Size> <centromere positions file: Chr \t centromere position> <genes gtf or bed> <output file>
#Author: Alice Naftaly, June 2020

import sys

#read chromosome sizes
#chromosome number is in roman numerals
#returns dictionary with key == chr number and value == size
def pull_chr_length():
    chr_length_file = sys.argv[1]
    chr_size_dict = {}
    with open(chr_length_file, 'r') as length_file:
        for line in length_file:
            if line.startswith("chr"):
                new_line = line.split()
                chr = new_line[0]
                size = int(new_line[1])
                chr_size_dict.update({chr:size})
    return chr_size_dict


#read centromere positions
#returns dictionary with key == chr number and value == centromere position
def pull_centromeres():
    centromeres_file = sys.argv[2]
    centromeres_dict = {}
    with open(centromeres_file, 'r') as centromeres:
        for line in centromeres:
            if line.startswith("chr"):
                new_line = line.split()
                chr = new_line[0]
                cent_position = int(new_line[1])
                centromeres_dict.update({chr:cent_position})
    return centromeres_dict


#create bins for chromosome
#arm 1 represents from the start of the chromosome to the centromere
#arm 2 represents from the centromere to the end of the chromosome
def create_bin_size():
    chr_length = pull_chr_length()
    centromere_pos = pull_centromeres()
    bins_dict = {}
    for chr in chr_length:
        if chr in centromere_pos:
            chr_size = chr_length[chr]
            centromere_position = centromere_pos[chr]
            arm_2 = chr_size - centromere_position
            arm_1 = chr_size - arm_2
            bin_size_arm1 = int(round(arm_1/10,0))
            bin_size_arm2 = int(round(arm_2/10,0))
            #bins variable will have [bin size arm 1, bin size arm 2]
            bin_sizes = [bin_size_arm1,bin_size_arm2]
            bins_dict.update({chr:bin_sizes})
    return bins_dict

#create bin boundaries
#returns dictionary with key == chr number and value == list of bins across chromosome
def create_bins_final():
    bin_sizes = create_bin_size()
    final_bins = {}
    for chr in bin_sizes:
        single_chr_bin = bin_sizes[chr]
        arm_1 = single_chr_bin[0]
        arm_2 = single_chr_bin[1]
        arm1_bin_count = 1
        arm1_position = 0
        while arm1_bin_count <= 10:
            arm_1_bin = [arm1_position, (arm_1*arm1_bin_count)-1]
            if chr in final_bins:
                final_bins[chr].append(arm_1_bin)
            elif chr not in final_bins:
                final_bins.update({chr:[arm_1_bin]})
            arm1_position += arm_1
            arm1_bin_count += 1
        arm2_bin_count = 1
        arm2_position = arm1_position
        while arm2_bin_count <= 10:
            arm_2_bin = [arm2_position, (arm1_position +  (arm_2*arm2_bin_count))-1]
            if chr in final_bins:
                final_bins[chr].append(arm_2_bin)
            elif chr not in final_bins:
                final_bins.update({chr:[arm_2_bin]})
            arm2_position += arm_2
            arm2_bin_count += 1
    return final_bins

#read gene or transcript positions
#need gtf or bed file for this
#if using a gtf file: will need to change feature
def pull_gene_positions_gtf():
    genes_file = sys.argv[3]
    genes_dict = {}
    with open(genes_file, 'r') as genes:
        for line in genes:
            new_line = line.split("\t")
            chr = new_line[0]
            feature = new_line[2]
            if feature == "gene":
                pos_1 = int(new_line[3])
                pos_2 = int(new_line[4])
                gene_info = new_line[8].split(" ")
                #will need to change the column used depending on the feature
                #gene = 1
                #transcript = 3
                gene_id = gene_info[1].strip(";")
                dict_value = [pos_1, pos_2, gene_id]
                if chr in genes_dict:
                    genes_dict[chr].append(dict_value)
                elif chr not in genes_dict:
                    genes_dict.update({chr:[dict_value]})
    return genes_dict

#need to write read bed file function and function to relate gene id with isoform id from isoseq



#splitting genes into bins based on positions
#if any genes overlap a bin boundary, will put these into a list to examine more closely to determine how they should be binned if need be
#there are 89 genes that do not fully fall within 1 bin (will just leave these out for now)
#returns dictionary with key == chr number and value == [index, [list of genes]]
def split_genes_into_bins():
    bins = create_bins_final()
    genes = pull_gene_positions_gtf()
    all_genes = []
    genes_in_bins = []
    genes_not_in_bins = []
    genes_in_bins_dict = {}
    for chr in bins:
        if chr in genes:
            single_chr_bins = bins[chr]
            single_chr_genes = genes[chr]
            for index,bin in enumerate(single_chr_bins):
                bin_start = bin[0]
                bin_end = bin[1]
                gene_list = []
                for single_gene in single_chr_genes:
                    gene_pos_1 = single_gene[0]
                    gene_pos_2 = single_gene[1]
                    gene_id = single_gene[2]
                    all_genes.append(gene_id)
                    if gene_pos_1 > bin_start and gene_pos_1 < bin_end and gene_pos_2 > bin_start and gene_pos_2 < bin_end:
                        gene_list.append(gene_id)
                        genes_in_bins.append(gene_id)
                dict_value = [index, gene_list]
                if chr in genes_in_bins_dict:
                    genes_in_bins_dict[chr].append(dict_value)
                elif chr not in genes_in_bins_dict:
                    genes_in_bins_dict.update({chr:[dict_value]})
    genes_not_in_bins = list(set(all_genes) - set(genes_in_bins))
    return genes_in_bins_dict


#write output file
def write():
    genes_in_bins = split_genes_into_bins()
    bins = create_bins_final()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for chr in genes_in_bins:
            single_chr_genes = genes_in_bins[chr]
            single_chr_bins = bins[chr]
            for index, value in enumerate(single_chr_genes):
                bin_start = single_chr_bins[index][0]
                bin_end = single_chr_bins[index][1]
                bin_number = value[0]
                bin_gene_list = value[1]
                for gene in bin_gene_list:
                    #chr number \t bin number \t bin start \t bin end \t gene id \n
                    final = "%s\t%s\t%s\t%s\t%s\n" % (str(chr), str(bin_number), str(bin_start), str(bin_end), str(gene))
                    out.write(final)


write()
