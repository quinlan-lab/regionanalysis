python addscores.py
python comparison.py > means.txt
awk 'FNR==NR{gerp=$1; phast=$2; cadd=$3} FNR!=NR{n=split($(NF-2),a,","); for (i in a) sum+=a[i]; if (sum/n < gerp && $12 > 90) print; sum=0}' means.txt $HOME/analysis/scoredregions.bed > lowgerp.txt
awk 'FNR==NR{gerp=$1; phast=$2; cadd=$3} FNR!=NR{n=split($(NF-1),a,","); for (i in a) sum+=a[i]; if (sum/n < phast && $12 > 90) print; sum=0}' means.txt $HOME/analysis/scoredregions.bed > lowphast.txt
awk 'FNR==NR{gerp=$1; phast=$2; cadd=$3} FNR!=NR{n=split($NF,a,","); for (i in a) sum+=a[i]; if (sum/n < cadd && $12 > 90) print; sum=0}' means.txt $HOME/analysis/scoredregions.bed > lowcadd.txt
