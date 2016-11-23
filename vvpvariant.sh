tr -s "" "\n" < ogfiles/VVP_KCNH2.txt | sed '1d' | awk '{printf $1"\t"$2-1"\t"} {for (i=2;i<=NF;i++) printf $i"\t"; printf "\n"}' > temp
#pathogenic variants in KCNH2
grep -v "Yes" temp > VVPvariants/pathogenicKCNH2.bed
rm temp
bedtools intersect -b <(sed '1d' exacresiduals/results/2016_11_17/weightedresiduals.txt) -a <(cut -f -24 VVPvariants/pathogenicKCNH2.bed) -wb | cut -f -24,36 | sort -k25,25n > VVPvariants/pathogenicKCNH2CCR.bed
echo "stats for pathogenic variants in KCNH2"
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/pathogenicKCNH2CCR.bed
#benign variants in KCNH2
echo "stats for benign clinvar variants in KCNH2"
cat <(grep "^#" benign.vcf) <(grep KCNH2 benign.vcf) > VVPvariants/benignKCNH2.vcf
bedtools intersect -a <(sed '1d' exacresiduals/results/2016_11_17/weightedresiduals.txt) -b VVPvariants/benignKCNH2.vcf | cut -f 12 > VVPvariants/benignKCNH2CCR.txt
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/benignKCNH2CCR.txt
#make coverage bed graph
awk '{print $1,$2-1,$2,$8}' exacresiduals/data/chr7.coverage.vcf | grep -v '^#' | awk '{split($4,a,","); print $1"\t"$2"\t"$3"\t"a[3]}' > chr7.cov.bedgraph
