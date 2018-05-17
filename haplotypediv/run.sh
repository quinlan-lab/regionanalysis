zcat $HOME/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz \
| cut -f 4,14 \
| sort \
| bedtools groupby -g 1 -c 2 -o max \
> ccr_gene_max.txt


###########################################################
#  https://academic.oup.com/hmg/article/26/3/489/2805387  #
###########################################################

#Non-Disease and Non-Essential
cat gene_essentiality.txt | awk '$2=="NDNE"' | cut -f1 > gene_essentiality.ndne.txt
#Essential Non-Disease
cat gene_essentiality.txt | awk '$2=="END"' | cut -f1 > gene_essentiality.end.txt
#Complex Non-Mendelian 
cat gene_essentiality.txt | awk '$2=="CNM"' | cut -f1 > gene_essentiality.cnm.txt
#Mendelian Non-Complex
cat gene_essentiality.txt | awk '$2=="MNC"' | cut -f1 > gene_essentiality.mnc.txt
#Complex-Mendelian
cat gene_essentiality.txt | awk '$2=="CM"' | cut -f1 > gene_essentiality.cm.txt

ISESS=$(python essentiality.py ccr_gene_max.txt gene_essentiality.end.txt | awk '$1>=99' | wc -l)
NOTESS=$(python essentiality.py ccr_gene_max.txt gene_essentiality.end.txt | awk '$1<99' | wc -l)
ISDIS=$(python essentiality.py ccr_gene_max.txt <( cat gene_essentiality.ndne.txt gene_essentiality.cm.txt gene_essentiality.cnm.txt gene_essentiality.mnc.txt) | awk '$1>=99' | wc -l)
NOTDIS=$(python essentiality.py ccr_gene_max.txt <( cat gene_essentiality.ndne.txt gene_essentiality.cm.txt gene_essentiality.cnm.txt gene_essentiality.mnc.txt) | awk '$1<99' | wc -l)

echo "Essentiality at 99% CCR"
python fisher.py $ISESS $NOTESS $ISDIS $NOTDIS

ISESS=$(python essentiality.py ccr_gene_max.txt gene_essentiality.end.txt | awk '$1>=95' | wc -l)
1003
NOTESS=$(python essentiality.py ccr_gene_max.txt gene_essentiality.end.txt | awk '$1<95' | wc -l)
510
ISDIS=$(python essentiality.py ccr_gene_max.txt <( cat gene_essentiality.ndne.txt gene_essentiality.cm.txt gene_essentiality.cnm.txt gene_essentiality.mnc.txt) | awk '$1>=95' | wc -l)
5681
NOTDIS=$(python essentiality.py ccr_gene_max.txt <( cat gene_essentiality.ndne.txt gene_essentiality.cm.txt gene_essentiality.cnm.txt gene_essentiality.mnc.txt) | awk '$1<95' | wc -l)
9367

echo "Essentiality at 95% CCR"
python fisher.py $ISESS $NOTESS $ISDIS $NOTDIS


###################################################################
#  https://link.springer.com/article/10.1007%2Fs00439-012-1243-6  #
###################################################################

ISSTRONG=$(python ld.py ccr_gene_max.txt strong_ld.txt | awk '$1>=99' | wc -l)
NOTSTRONG=$(python ld.py ccr_gene_max.txt strong_ld.txt | awk '$1<99' | wc -l)
ISWEAK=$(python ld.py ccr_gene_max.txt weak_ld.txt | awk '$1>=99' | wc -l)
NOTWEAK=$(python ld.py ccr_gene_max.txt weak_ld.txt | awk '$1<99' | wc -l)

echo "LD at 99%"
python fisher.py $ISSTRONG $NOTSTRONG $ISWEAK $NOTWEAK

ISSTRONG=$(python ld.py ccr_gene_max.txt strong_ld.txt | awk '$1>=95' | wc -l)
NOTSTRONG=$(python ld.py ccr_gene_max.txt strong_ld.txt | awk '$1<95' | wc -l)
ISWEAK=$(python ld.py ccr_gene_max.txt weak_ld.txt | awk '$1>=95' | wc -l)
NOTWEAK=$(python ld.py ccr_gene_max.txt weak_ld.txt | awk '$1<95' | wc -l)

echo "LD at 99%"
python fisher.py $ISSTRONG $NOTSTRONG $ISWEAK $NOTWEAK
