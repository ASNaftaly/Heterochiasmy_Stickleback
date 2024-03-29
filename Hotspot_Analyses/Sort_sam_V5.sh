ml SAMtools/1.9-foss-2016b

#samtools view -bS LW_hotspots_V5_positions.sam > LW_hotspots_V5_positions.bam
#samtools sort LW_hotspots_V5_positions.bam > LW_hotspots_V5_positions_sorted.bam
#samtools index LW_hotspots_V5_positions_sorted.bam

#samtools view -bS PS_hotspots_V5_positions.sam > PS_hotspots_V5_positions.bam
#samtools sort PS_hotspots_V5_positions.bam > PS_hotspots_V5_positions_sorted.bam
#samtools index PS_hotspots_V5_positions_sorted.bam

samtools view -bS LW_shared_hotspots_V5_positions.sam > LW_shared_hotspots_V5_positions.bam
samtools sort LW_shared_hotspots_V5_positions.bam > LW_shared_hotspots_V5_positions_sorted.bam
samtools index LW_shared_hotspots_V5_positions_sorted.bam

samtools view -bS PS_shared_hotspots_V5_positions.sam > PS_shared_hotspots_V5_positions.bam
samtools sort PS_shared_hotspots_V5_positions.bam > PS_shared_hotspots_V5_positions_sorted.bam
samtools index PS_shared_hotspots_V5_positions_sorted.bam

samtools view -bS LW_unique_hotspots_V5_positions.sam > LW_unique_hotspots_V5_positions.bam
samtools sort LW_unique_hotspots_V5_positions.bam > LW_unique_hotspots_V5_positions_sorted.bam
samtools index LW_unique_hotspots_V5_positions_sorted.bam

samtools view -bS PS_unique_hotspots_V5_positions.sam > PS_unique_hotspots_V5_positions.bam
samtools sort PS_unique_hotspots_V5_positions.bam > PS_unique_hotspots_V5_positions_sorted.bam
samtools index PS_unique_hotspots_V5_positions_sorted.bam
