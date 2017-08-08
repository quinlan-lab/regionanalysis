# using ad_genes.bed from pathoscore's make.sh for AD genes
# AD genes covered by top 10%, 5% and 1% CCRs
bedtools intersect -a ad_genes.bed -b <(awk '$14>=90' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a ad_genes.bed -b <(awk '$14>=95' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a ad_genes.bed -b <(awk '$14>=99' <(gzip -dc gnomad-ccrs.bed.gz)) | awk '{t+=$3-$2} END {print t}'
# total number of bp in ad_genes
awk '{t+=$3-$2} END {print t}' ad_genes.bed
#bp of ad_genes covered by %X CCR/total bp of AD genes gives denominator
#number of pathogenic variants in AD genes covered by top 10%, 5% and 1% CCRs
bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=90' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l
bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=95' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l
bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b <(awk '$14>=99' <(gzip -dc gnomad-ccrs.bed.gz)) -u -header | bedtools intersect -a stdin -b ad_genes.bed -u | wc -l
#total number of pathogenic variants in AD genes
bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -b ad_genes.bed -u | wc -l
#number of AD pathogenic variants covered by %X CCR/total number of AD pathogenic variants gives numerator
#divide numerator by denominator to get result for each % CCR
