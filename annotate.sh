#download GENE ID mapped to Associated Gene Name from Grch37 biomart == mart_export.txt
cp exacresiduals/results/2017_04_09/weightedresiduals-cpg-novariant.txt weightedresiduals.bed
awk 'FNR==NR{genes[$2]=$1;next} {for (i in genes) if (i==$4) print $0"\t"genes[i]}' mart_export.txt weightedresiduals.bed | sort -k1,1 -k2,2n | bgzip -c > ccrs.bed.gz
tabix ccrs.bed.gz
#ccr.conf refers to appropriate columns and file
vcfanno -lua $HOME/software/vcfanno/docs/examples/clinvar_exac.lua -p 15 -base-path . ccr.conf $DATA/clinvar-vt-anno-vep.vcf.gz > clinvar_ccr.vcf
cat <(grep '^#' clinvar_ccr.vcf) <(grep -v '^#' clinvar_ccr.vcf | sort -k1,1 -k2,2n) | bgzip -c > clinvar_ccr.vcf.gz
tabix clinvar_ccr.vcf.gz
