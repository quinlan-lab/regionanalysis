# http://www.cell.com/cms/attachment/2113541342/2084322352/mmc1.pdf (Table S14) ignoring the X-linked dominant, 11 genes with mutational clustering in DD and ID patients and no HI effects (from AJHG: http://www.cell.com/ajhg/fulltext/S0002-9297(17)30326-9) #used DECIPHER (DDG2P) webtool for clustering
#table S14
awk 'FNR==NR{genes[$1]; next} {for (gene in genes) if (gene == $4) print}' ogfiles/nhigenes.txt <(zcat $HOME/analysis/essentials/gnomadbased-ccrs.bed.gz) | awk '$NF >= 99' > nhiddregions.bed
echo "Number of unique regions >=99 in 11 autosomal NHI DD genes with mutational clusters"
cut -f 4,14 nhiddregions.bed | sort | uniq | wc -l
echo "Number of genes covered by regions >=99 in 11 autosomal NHI DD genes with mutational clusters"
cut -f 4 nhiddregions.bed | sort | uniq | wc -l

#table 1
awk 'FNR==NR{genes[$1]; next} {for (gene in genes) if (gene == $4) print}' ogfiles/ddgenes.txt <(zcat $HOME/analysis/essentials/gnomadbased-ccrs.bed.gz) | awk '$NF >= 99' > ddregions.bed
echo "Number of unique regions >=99 in 14 autosomal DD genes with mutational clusters"
cut -f 4,14 ddregions.bed | sort | uniq | wc -l
echo "Number of genes covered by regions >=99 in 14 autosomal DD genes with mutational clusters"
cut -f 4 ddregions.bed | sort | uniq | wc -l
