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

if [ ! -s $DATA/ExAC.r1.vt.vep.alts.vcf.gz ]; then
    cat <(grep "^#" $DATA/ExAC.r1.vt.vep.vcf) <(grep -v "^#" $DATA/ExAC.r1.vt.vep.vcf | awk '{print $0 ";refs=" $4 ";" "alts=" $5}' FS='\t' OFS='\t') > $DATA/ExAC.r1.vt.vep.alts.vcf; bgzip -c $DATA/ExAC.r1.vt.vep.alts.vcf > $DATA/ExAC.r1.vt.vep.alts.vcf.gz; tabix $DATA/ExAC.r1.vt.vep.alts.vcf.gz
fi

if [ ! -s $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.alts.vcf]; then
    cat <(grep "^#" $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.vcf) <(grep -v "^#" $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.vcf | awk '{print $0 ";refs=" $4 ";" "alts=" $5}' FS='\t' OFS='\t') > $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.alts.vcf; bgzip -c $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.alts.vcf > $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.alts.vcf.gz; tabix $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.alts.vcf.gz
fi

# cat <(grep "^#" $DATA/gnomad_vt.vcf) <(grep -v "^#" $DATA/gnomad_vt.vcf | awk '{print $0 ";refs=" $4 ";" "alts=" $5}' FS='\t' OFS='\t') command to add refs and alts here, then add gnomad_alts from alts and exac_alts from alts to conf files depending on file annotated with

if [ $exac ]; then
    vcfanno -lua $HOME/analysis/custom.lua -permissive-overlap -p 4 -base-path $DATA $HOME/analysis/exaconly.conf $DATA/${FILE}_vt.vcf > $DATA/${FILE}-anno-vt.vcf
else
    vcfanno -lua $HOME/analysis/custom.lua -permissive-overlap -p 4 -base-path $DATA $HOME/analysis/annovars.conf $DATA/${FILE}_vt.vcf > $DATA/${FILE}-anno-vt.vcf
fi

perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/${FILE}-anno-vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/${FILE}-vep-anno-vt.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite

bgzip -c $DATA/${FILE}-vep-anno-vt.vcf > $DATA/${FILE}-vep-anno-vt.vcf.gz; tabix $DATA/${FILE}-vep-anno-vt.vcf.gz
