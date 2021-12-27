ml Python

for ((n=0;n<10000;n++))
do
  python3 Random_expectation_peak_distribution.py <peaks file> Chromosome_Sizes.txt <output random peak file>
  python3 Make20Bins_random_atac.py CentromerePositions_roman.txt <output random peak file> Chromosome_Sizes.txt <final randomized peaks file; will be appended to for 10k permutations>
  rm <output random peak file>
  echo $n
done
