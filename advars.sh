echo "CCRs 95 and above harboring at least one pathogenic variant"
awk '$NF>=95' patho-ccr.txt | grep -v "_exclude" | cut -f 9,12,15 | sort -k1,1 -k2,2 -k3,3 | uniq | wc -l
zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=95' | cut -f 1,4,7 | sort -k1,1 -k2,2 -k3,3 | uniq | wc -l
