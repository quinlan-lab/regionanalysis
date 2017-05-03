exac=$2 # if exaconly, put something here
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

if [ $exac ]; then
    vcfanno -lua $HOME/analysis/custom.lua -permissive-overlap -p 4 -base-path $DATA $HOME/analysis/exaconly.conf $DATA/${FILE}_vt.vcf > $DATA/${FILE}-anno-vt.vcf
else
    vcfanno -lua $HOME/analysis/custom.lua -permissive-overlap -p 4 -base-path $DATA $HOME/analysis/annovars.conf $DATA/${FILE}_vt.vcf > $DATA/${FILE}-anno-vt.vcf
fi

perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/${FILE}-anno-vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/${FILE}-vep-anno-vt.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite

bgzip -c $DATA/${FILE}-vep-anno-vt.vcf > $DATA/${FILE}-vep-anno-vt.vcf.gz; tabix $DATA/${FILE}-vep-anno-vt.vcf.gz
