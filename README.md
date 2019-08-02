# Heterochiasmy_Stickleback
Scripts Associated with Heterochiasmy Analyses in Threespine Stickleback Fish

# Project Design

I did ATAC-seq and RNA-seq on paired samples for two different populations of threespine stickleback fish (Pacific Ocean, a marine population, and Lake Washington, a freshwater population) in both males and females. I collected liver and gonads (ovaries or testes) from each sample where half the liver and one gonad were used for ATAC-seq and the other half of the liver and the other gonad were used for RNA-seq.

The scripts for calling peaks for ATAC-seq samples and the general pipeline for RNA-seq alignment to differential expression are in separate repositories (maybe, I might change my mind on this).


# After calling ATAC-seq Peaks
### 1. Getting Read Coverage per bp:
   GenomeCoverage.sh --> CalcChrLength.sh
   
### 2. Splitting Peaks into 20 bins per chromosome (to determine if peaks are uniformly distributed across chromosomes or if there are any tissue specific patterns):
   Make20Bins_ATAC_peaks.py == observed bin values
   (Random_expectation_peak_distribution.py --> Make20Bins_random_atac.py)x10k --> Random_peak_dist_summary.py -->        CalculateP_randompeakdist.py == expected peak distribution which is then compared to observed bin values to calculate significance for each bin
   
### 3. Gene Enrichment (Are peaks enriched around genes?):
   Peak_tssenrichment_500bp.py == observed peaks within 500bp of TSS
   (Tss_expected_enrich_ATACpeaks_500bp.py)x10k == creates random distribution to compare observed peaks near TSSs to and calculate significance
### 4. GO term
   PullGeneIDs_atac.py (optional --> CompareIDS_atac.py == compare samples to get unique/shared gene lists) --> Pull GO terms from Biomart --> PullGOterms_atac.py --> GoTermCount_atac.py == observed GO terms
   (Random_GO_distribution_atac.py)x10k --> Summarize_GO_permutations_atac.py --> SigGOterms_atac.py == creates random distribution to compare observed GO terms to and calcualtes significance
