#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=4:00:00
python residualgist/exac-resid.py > $DATA/exac-resid.txt
#python residualgist/resid-plot.py $DATA/exac-resid.txt > residualgist/foo.txt
#cat <(head -1 residualgist/foo.txt) <(sed '1d' residualgist/foo.txt | bash dups.sh) > residualgist/residuals.txt; rm residualgist/foo.txt
FILE=$1;FOLDER=$2
#sed '1d' residualgist/residuals.txt | awk '$9 >= 99' > regions/topresid.txt
#python middle.py > regions/midresid.txt
#awk 'NR==FNR{l+=$3-$2;next} {s+=$3-$2; if (s>=l) exit; else print}' regions/topresid.txt <(sed '1d' residualgist/residuals.txt | awk '$9 <= 10' | sort -k9,9n) > regions/bottomresid.txt
#sed '1d' residualgist/residuals.txt | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3,8 -o distinct,min,max,mean | awk '{print $2,$3,$4,$1,$5}' OFS="\t" | awk '{gene[$4]=$5; row[$4]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],100*rank[i]/NR}' > meangeneresiduals.txt
#sed '1d' residualgist/residuals.txt | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3,8 -o distinct,min,max,max | awk '{print $2,$3,$4,$1,$5}' OFS="\t" | awk '{gene[$4]=$5; row[$4]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],100*rank[i]/NR}' > maxgeneresiduals.txt
#awk '$6 > 99' meangeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > regions/topmeangeneresid.txt
#awk '$5 > -0.033 && $5 < 0.033' meangeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > regions/midmeangeneresid.txt
#awk '$6 > 99' maxgeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > regions/topmaxgeneresid.txt
#awk '$6 < 5' maxgeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > regions/midmaxgeneresid.txt
#bash loop.sh $FILE
#bash plotloop.sh $FOLDER
#bash analyze.sh 0 regions/bottomgerpresid.txt ogfiles/all_ar.tsv
#python plot.py results/all_ar.overlap.txt
