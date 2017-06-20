original=$1
fileName="${original##*/}"
fileExt=${fileName#*.}
FILE=${fileName%*.$fileExt}
if [ ! -s $DATA/${FILE}_dec.vcf ]; then
    vt decompose $original -o $DATA/${FILE}_dec.vcf -s 
fi

if [ ! -s $DATA/${FILE}_vt.vcf ]; then
    vt normalize $DATA/${FILE}_dec.vcf -o $DATA/${FILE}_vt.vcf -r $DATA/grch37.fa
fi

perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/${FILE}_vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/${FILE}-vep-vt.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite

sort -k1,1 -k2,2n $DATA/${FILE}-vep-vt.vcf | bgzip -c > $DATA/${FILE}-vep-vt.vcf.gz; tabix $DATA/${FILE}-vep-vt.vcf.gz
