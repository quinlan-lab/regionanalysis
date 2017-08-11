# using ad_genes.bed from pathoscore's make.sh for AD genes
# AD genes covered by top 10%, 5% and 1% CCRs
## GET GENES FROM THIS URL and store in genescreens/adstrict_genecards.txt
#http://www.genecards.org/Search/Keyword?queryString=%20%5Bdisorders%5D%20%20(%20%20autosomal%20AND%20dominant%20%20)%20%20AND%20%20%5Bvariants%5D%20%20(%20%20autosomal%20AND%20dominant%20%20)%20&pageSize=-1&startPage=0

#number of AD pathogenic variants covered by %X CCR/total number of AD pathogenic variants gives numerator
#bp of ad_genes covered by %X CCR/total bp of AD genes gives denominator
#divide numerator by denominator to get result for each % CCR

# total number of bp in ad_genes
TBP=$(awk '{t+=$3-$2} END {print t}' ad_genes.bed)

#total number of pathogenic variants in AD genes
TPV=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b ad_genes.bed -u | wc -l)

BP=$(bedtools intersect -a ad_genes.bed -b <(awk '$14>=90' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}') # bp covered in AD genes by >=90% CCRs
PV=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=90' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l) # vars covered in AD genes by >=90% CCRs
echo "90% CCR pathogenic enrichment in AD genes"
echo "scale=5; ($PV/$TPV)/($BP/$TBP)" | bc
BP=$(bedtools intersect -a ad_genes.bed -b <(awk '$14>=95' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}')
PV=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=95' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l)
echo "95% CCR pathogenic enrichment in AD genes"
echo "scale=5; ($PV/$TPV)/($BP/$TBP)" | bc
BP=$(bedtools intersect -a ad_genes.bed -b <(awk '$14>=99' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}')
PV=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=99' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l)
echo "99% CCR pathogenic enrichment in AD genes"
echo "scale=5; ($PV/$TPV)/($BP/$TBP)" | bc
BP=$(bedtools intersect -a ad_genes.bed -b <(awk '$14>=99.5' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}')
PV=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=99.5' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l)
echo "99.5% CCR pathogenic enrichment in AD genes"
echo "scale=5; ($PV/$TPV)/($BP/$TBP)" | bc
BP=$(bedtools intersect -a ad_genes.bed -b <(awk '$14>=99.9' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}')
PV=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=99.9' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l)
echo "99.9% CCR pathogenic enrichment in AD genes"
echo "scale=5; ($PV/$TPV)/($BP/$TBP)" | bc
