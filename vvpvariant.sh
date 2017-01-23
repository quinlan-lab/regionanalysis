tr -s "" "\n" < ogfiles/VVP_KCNH2.txt | sed '1d' | awk '{printf $1"\t"$2-1"\t"} {for (i=2;i<=NF;i++) printf $i"\t"; printf "\n"}' > /tmp/temp
#pathogenic variants in KCNH2
grep -v "Yes" /tmp/temp > VVPvariants/pathogenicKCNH2.bed
grep -v "No" /tmp/temp > VVPvariants/benignKCNH2.bed
bedtools intersect -a <(cut -f -24 VVPvariants/pathogenicKCNH2.bed) -b <(sed '1d' exacresiduals/results/2016_12_01/weightedresiduals.txt) -wb | cut -f -24,36 | sort -k25,25n > VVPvariants/pathogenicKCNH2CCR.bed
bedtools intersect -a <(cut -f -24 VVPvariants/benignKCNH2.bed) -b <(sed '1d' exacresiduals/results/2016_12_01/weightedresiduals.txt) -wb | cut -f -24,36 | sort -k25,25n > VVPvariants/benignKCNH2CCR.bed
echo "stats for pathogenic zebrafish variants in KCNH2"
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/pathogenicKCNH2CCR.bed
#benign variants in KCNH2
echo "stats for benign zebrafish variants in KCNH2"
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/benignKCNH2CCR.bed
echo "stats for benign clinvar variants in KCNH2"
cat <(grep "^#" benign.vcf) <(grep KCNH2 benign.vcf) > VVPvariants/benignKCNH2.vcf
bedtools intersect -a <(sed '1d' exacresiduals/results/2016_12_01/weightedresiduals.txt) -b VVPvariants/benignKCNH2.vcf | cut -f 12 > VVPvariants/benignKCNH2CCR.txt
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/benignKCNH2CCR.txt
#pathogenic variants in LQTS
grep "pathogenic" VVPvariants/LQTSvariants.bed | sed 's/^chr//g' > VVPvariants/pathogenicLQTS.bed
bedtools intersect -a VVPvariants/pathogenicLQTS.bed -b <(sed '1d' exacresiduals/results/2016_12_01/weightedresiduals.txt) -wb | cut -f -7,19 | sort -k8,8n > VVPvariants/pathogenicLQTSCCR.bed
echo "stats for pathogenic variants in LQTS"
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/pathogenicLQTSCCR.bed
#benign variants in LQTS
grep "Benign" VVPvariants/LQTSvariants.bed | sed 's/^chr//g' > VVPvariants/benignLQTS.bed
bedtools intersect -a VVPvariants/benignLQTS.bed -b <(sed '1d' exacresiduals/results/2016_12_01/weightedresiduals.txt) -wb | cut -f -7,19 | sort -k8,8n > VVPvariants/benignLQTSCCR.bed
echo "stats for benign variants in LQTS"
awk 'NR == 1{ max=$NF; min=$NF; sum=0 }
     {arr[NR]=$NF} { if ($NF>max) max=$NF; if ($NF<min) min=$NF; sum+=$NF}
     END {printf "Min: %d\tMax: %d\tAverage: %f", min, max, sum/NR}
     END { if (NR%2==1) printf " Median: " arr[(NR+1)/2] "\n";
     else print " Median: " (arr[NR/2]+arr[NR/2+1])/2 "\n"}' VVPvariants/benignLQTSCCR.bed
#make coverage bed graph
awk '{print $1,$2-1,$2,$8}' exacresiduals/data/chr7.coverage.vcf | grep -v '^#' | awk '{split($4,a,","); print $1"\t"$2"\t"$3"\t"a[3]}' > chr7.cov.bedgraph
