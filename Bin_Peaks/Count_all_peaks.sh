#shell script to count total number of peaks for xls file

for file in O1_filteredpeaks.xls
do
  export sample=${file%%peaks.xls}
  echo $sample
  chr1=$(grep -i -P "chrI\t" $file | wc -l)
  chr2=$(grep -i -P "chrII\t" $file | wc -l)
  chr3=$(grep -i -P "chrIII\t" $file | wc -l)
  chr4=$(grep -i -P "chrIV\t" $file | wc -l)
  chr5=$(grep -i -P "chrV\t" $file | wc -l)
  chr6=$(grep -i -P "chrVI\t" $file | wc -l)
  chr7=$(grep -i -P "chrVII\t" $file | wc -l)
  chr8=$(grep -i -P "chrVIII\t" $file | wc -l)
  chr9=$(grep -i -P "chrIX\t" $file | wc -l)
  chr10=$(grep -i -P "chrX\t" $file | wc -l)
  chr11=$(grep -i -P "chrXI\t" $file | wc -l)
  chr12=$(grep -i -P "chrXII\t" $file | wc -l)
  chr13=$(grep -i -P "chrXIII\t" $file | wc -l)
  chr14=$(grep -i -P "chrXIV\t" $file | wc -l)
  chr15=$(grep -i -P "chrXV\t" $file | wc -l)
  chr16=$(grep -i -P "chrXVI\t" $file | wc -l)
  chr17=$(grep -i -P "chrXVII\t" $file | wc -l)
  chr18=$(grep -i -P "chrXVIII\t" $file | wc -l)
  chr19=$(grep -i -P "chrXIX\t" $file | wc -l)
  chr20=$(grep -i -P "chrXX\t" $file | wc -l)
  chr21=$(grep -i -P "chrXXI\t" $file | wc -l)
  chrY=$(grep -i -P "chrY\t" $file | wc -l)
  echo $chr1
  echo $chr2
  echo $chr3
  echo $chr4
  echo $chr5
  echo $chr6
  echo $chr7
  echo $chr8
  echo $chr9
  echo $chr10
  echo $chr11
  echo $chr12
  echo $chr13
  echo $chr14
  echo $chr15
  echo $chr16
  echo $chr17
  echo $chr18
  echo $chr20
  echo $chr21
  echo $chr19
  echo $chrY
  autosome_total=$(($chr1 + $chr2 + $chr3 + $chr4 + $chr5 + $chr6 + $chr7 + $chr8 + $chr9 + $chr10 + $chr11 + $chr12 + $chr13 + $chr14 + $chr15 + $chr16 + $chr17 + $chr18 + $chr20 + $chr21))
  final_total_count=$(($chr1 + $chr2 + $chr3 + $chr4 + $chr5 + $chr6 + $chr7 + $chr8 + $chr9 + $chr10 + $chr11 + $chr12 + $chr13 + $chr14 + $chr15 + $chr16 + $chr17 + $chr18 + $chr19 + $chr20 + $chr21 + $chrY))
  echo $autosome_total
  echo $final_total_count
done
