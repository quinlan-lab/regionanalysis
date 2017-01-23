vt decompose $DATA/clinvar_20170104.vcf.gz -o $DATA/clinvar_20170104_dec.vcf

vt normalize $DATA/clinvar_20170104_dec.vcf -o $DATA/clinvar_vt.vcf -r $DATA/grch37.fa

vcfanno -lua $HOME/software/vcfanno/docs/examples/clinvar_exac.lua -p 4 -base-path $DATA $HOME/software/vcfanno/docs/examples/clinvar_exac.conf $DATA/clinvar_vt.vcf > $DATA/clinvar-anno-vt.vcf

perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/clinvar-anno-vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/clinvar-vt-anno-vep.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM --offline

bgzip $DATA/clinvar-vt-anno-vep.vcf; tabix $DATA/clinvar-vt-anno-vep.vcf.gz
