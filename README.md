# Heterochiasmy_Stickleback
<<<<<<< HEAD
#Scripts Associated with Heterochiasmy analyses in Threespine Stickleback
=======
Scripts Associated with Heterochiasmy Analyses in Threespine Stickleback Fish

# Project Design

I did ATAC-seq and RNA-seq on paired samples for two different populations of threespine stickleback fish (Pacific Ocean, a marine population, and Lake Washington, a freshwater population) in both males and females. I collected liver and gonads (ovaries or testes) from each sample where half the liver and one gonad were used for ATAC-seq and the other half of the liver and the other gonad were used for RNA-seq.

The scripts for calling peaks for ATAC-seq samples and the general pipeline for RNA-seq alignment to differential expression are in separate repositories (maybe, I might change my mind on this).


# After calling ATAC-seq Peaks
These analyses can be run in any order

### Getting Read Coverage per bp:
   GenomeCoverage.sh
  
### Creating Chromosome Length File
   CalcChrLength.sh
   
### Splitting Peaks into 20 bins per chromosome (to determine if peaks are uniformly distributed across chromosomes or if there are any tissue specific patterns):
   Make20Bins_ATAC_peaks.py == observed bin values
   
   (Random_expectation_peak_distribution.py --> Make20Bins_random_atac.py)x10k --> Random_peak_dist_summary.py -->        CalculateP_randompeakdist.py == expected peak distribution which is then compared to observed bin values to calculate significance for each bin
   
### Gene Enrichment (Are peaks enriched around genes?):
   Peak_tssenrichment_500bp.py == observed peaks within 500bp of TSS
   
   (Tss_expected_enrich_ATACpeaks_500bp.py)x10k == creates random distribution to compare observed peaks near TSSs to and calculate significance
### GO term enrichment
   PullGeneIDs_atac.py (optional --> CompareIDS_atac.py == compare samples to get unique/shared gene lists) --> Pull GO terms from Biomart --> PullGOterms_atac.py --> GoTermCount_atac.py == observed GO terms
   
   (Random_GO_distribution_atac.py)x10k --> Summarize_GO_permutations_atac.py --> SigGOterms_atac.py == creates random distribution to compare observed GO terms to and calcualtes significance

### Transposable Element enrichment
   *I separated these TEs to suit my needs. I broadly examined TEs as class 1 = RNA(LINES,SINES, LTRs), class II = DNA, and RC/helitrons (I've seen these under class II TEs or even in their own class, so I examined them separately). I also examined TEs more specifically by looking at TE families under each broad class*
   
   Examining TEs separated as DNA, LINES, SINES, LTRs, RC/Helitrons
   TE_observed_counts_ATACpeaks.py == observed TE counts for 5 large TE classes
   
   (<BroadTEclass>TE_expected_counts_ATAC_peaks.py)x10k --> PullSum_TE_expect.py == randomizes peaks to determine which TE classes are significantly enriched
   
   Examining TEs at a family level
   Pull_TEfamilies_at_peaks.py --> TEfamily_counts.py == counts observed peaks at each TE family
   
   (TEfamily_random_expectations.py --> TEfamily_counts.py)x10k --> Summarize_TEfamily_permutations.py --> CalculateP_randomTEfamily.py == ranomdizes peaks to determine which TE familes are significantly enriched
>>>>>>> 518218f2c87fb0e3483f6f4efe21161550c5da2c
