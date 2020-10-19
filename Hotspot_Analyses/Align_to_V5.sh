module load Bowtie2/2.3.5.1-foss-2018a
export btBuild=V5_Genome_Assembly/Bowtie2_index/StickV5_bowtie2_index

bowtie2 -p 8 -x $btBuild -f LW_hotspots_glazer_positions.fasta -S LW_hotspots_V5_positions.sam

bowtie2 -p 8 -x $btBuild -f PS_hotspots_glazer_positions.fasta -S PS_hotspots_V5_positions.sam

bowtie2 -p 8 -x $btBuild -f LW_shared_hotspots_glazer_positions.fasta -S LW_shared_hotspots_V5_positions.sam

bowtie2 -p 8 -x $btBuild -f PS_shared_hotspots_glazer_positions.fasta -S PS_shared_hotspots_V5_positions.sam

bowtie2 -p 8 -x $btBuild -f LW_unique_hotspots_glazer_positions.fasta -S LW_unique_hotspots_V5_positions.sam

bowtie2 -p 8 -x $btBuild -f PS_unique_hotspots_glazer_positions.fasta -S PS_unique_hotspots_V5_positions.sam
