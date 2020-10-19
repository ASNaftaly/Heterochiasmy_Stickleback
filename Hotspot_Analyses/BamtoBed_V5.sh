#PBS -S /bin/bash
#PBS -q batch
#PBS -N convert_bam_to_bed_V5
#PBS -l nodes=1:ppn=1
#PBS -l mem=15gb
#PBS -l walltime=:00:00
#PBS -M alice.shanfelter@uga.edu
#PBS -m abe

cd /scratch/afs16076/Paired_ATAC_RNA/ATAC_seq/Hotspot_comparisons/

module load BEDTools/2.29.2-GCC-8.2.0-2.31.1

#bedtools bamtobed -i LW_hotspots_V5_positions_sorted.bam >> LW_hotspots_V5_positions.bed

#bedtools bamtobed -i PS_hotspots_V5_positions_sorted.bam >> PS_hotspots_V5_positions.bed

bedtools bamtobed -i LW_shared_hotspots_V5_positions.bam >> LW_shared_hotspots_V5_positions.bed

bedtools bamtobed -i PS_shared_hotspots_V5_positions.bam >> PS_shared_hotspots_V5_positions.bed

bedtools bamtobed -i LW_unique_hotspots_V5_positions.bam >> LW_unique_hotspots_V5_positions.bed

bedtools bamtobed -i PS_unique_hotspots_V5_positions.bam >> PS_unique_hotspots_V5_positions.bed
