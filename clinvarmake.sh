exacver=$1 # conf file for vcfanno, either clinvar_gnomad.conf or clinvar_exac.conf

if [ ! -s $DATA/clinvar_20170104_dec.vcf ]; then
    vt decompose $DATA/clinvar_20170104.vcf.gz -o $DATA/clinvar_20170104_dec.vcf
fi

if [ ! -s $DATA/clinvar_vt.vcf ]; then
    vt normalize $DATA/clinvar_20170104_dec.vcf -o $DATA/clinvar_vt.vcf -r $DATA/grch37.fa
fi

vcfanno -lua custom.lua -p 4 -base-path $DATA clinvar_gnomad.conf $DATA/clinvar_vt.vcf > $DATA/clinvar-anno-vt.vcf

perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/clinvar-anno-vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/clinvar-vt-anno-vep.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM --offline --fork 12 --force_overwrite

bgzip -c $DATA/clinvar-vt-anno-vep.vcf > $DATA/clinvar-vt-anno-vep.vcf.gz; tabix $DATA/clinvar-vt-anno-vep.vcf.gz
