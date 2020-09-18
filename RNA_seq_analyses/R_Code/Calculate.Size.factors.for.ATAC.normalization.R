library(DESeq2)
library(tximport)

#using all genes

#transcript to gene file
transcript2gene <- read.table("Isoform.IDs.to.Gene.IDs.txt", header=T)

#reading in kallisto files for all genes
files <- file.path("/Users/Alice/Documents/UGA/White_Lab/Project/Paired_ATAC_RNA/RNA_seq/Expression_Analyses", c("l1f_rna_abundance.tsv", "l2f_rna_abundance.tsv", "l1m_rna_abundance.tsv", "l2m_rna_abundance.tsv", "o1_rna_abundance.tsv", "o2_rna_abundance.tsv", "t1_rna_abundance.tsv","t2_rna_abundance.tsv"))
names(files) <- c("L1F", "L2F", "L1M", "L2M", "O1", "O2", "T1", "T2")
all_txi <- tximport(files, type="kallisto",tx2gene=transcript2gene)

samples <- read.table("Samples_info.txt",header=T)
rownames(samples) <- samples$Sample
samples[,c("Sample","Tissue")]

#create dds object for all genes
dds <- DESeqDataSetFromTximport(all_txi,colData=samples,design=~Tissue)
dds <- estimateSizeFactors(dds)
#this gives the normalized factors 
sf  <- estimateSizeFactorsForMatrix(normalizationFactors(dds)) 



