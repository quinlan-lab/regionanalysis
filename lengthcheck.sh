sed '1d' exacresiduals/fullresiduals.txt | cut -f 5 | sort | uniq > /tmp/transcripts.txt
cut -f 5 regions/topresid.txt | sort | uniq > /tmp/toptranscripts.txt
cut -f 4,5,9 /scratch/ucgd/lustre/u1021864/serial/vepcanonicalexons.gtf | cut -d " " -f 1,4 | perl -pe 's/gene_id|\"|\;//g' | awk '{l=$2-$1+1} {print $3,l}' | tr -s ' ' '\t' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum > /tmp/transcriptlengths.txt
awk 'NR==FNR{a[$1]} $1 in a {t+=$2} END {print t}' /tmp/transcripts.txt /tmp/transcriptlengths.txt
