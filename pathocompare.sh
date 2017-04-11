#for generating ladder plots and ROC curves
#first run bash clinvarmake.sh if the clinvar file below is not yet generated
TITLE=$1

python clinvarfilter.py -c $DATA/clinvar-vt-anno-vep.vcf.gz -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria

cat <(grep '^#' benign.vcf) <(grep -v '^#' benign.vcf | sort -k1,1 -k2,2n) | bgzip -c > benign.vcf.gz; tabix benign.vcf.gz
cat <(grep '^#' patho.vcf) <(grep -v '^#' patho.vcf | sort -k1,1 -k2,2n) | bgzip -c > patho.vcf.gz; tabix patho.vcf.gz

#exac
bedtools intersect -a <(sed '1d' exacresiduals/results/2017_03_30/weightedresiduals-cpg-novariant.txt) -b patho.vcf.gz | cut -f 14 > tmp/ccrpatho
bedtools intersect -a <(sed '1d' exacresiduals/results/2017_03_30/weightedresiduals-cpg-novariant.txt) -b benign.vcf.gz | cut -f 14 > tmp/ccrbenign

#gnomAD
bedtools intersect -a <(sed '1d' exacresiduals/results/2017_04_09/weightedresiduals-cpg-novariant.txt) -b patho.vcf.gz | cut -f 14 > tmp/ccr2patho
bedtools intersect -a <(sed '1d' exacresiduals/results/2017_04_09/weightedresiduals-cpg-novariant.txt) -b benign.vcf.gz | cut -f 14 > tmp/ccr2benign

# generate bed file from pli file (and a vcf CADD file by manually adding a header to the 1.3 version TSV file)

if [ ! -s pli.bed ]; then
    bedtools groupby -i $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt -g 3,5,6,2 -c 20 -o collapse | tr -s ' ' '\t' > pli.bed
fi
if [ ! -s $DATA/CADD.vcf.gz ]; then
    cat caddheader <(sed '1,2d' $DATA/whole_genome_SNVs.tsv | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS="\t") > $DATA/CADD.vcf # whole_genome_SNVs.tsv is the original CADD file
    bgzip -@ 12 -c $DATA/CADD.vcf > $DATA/CADD.vcf.gz; tabix $DATA/CADD.vcf.gz
fi

if [ ! -s $DATA/CADDindels.vcf.gz ]; then
    cat caddheader <(sed '1,2d' $DATA/InDels.tsv | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS="\t") > $DATA/CADDindels.vcf
    bgzip -@ 12 -c $DATA/CADDindels.vcf > $DATA/CADDindels.vcf.gz; tabix $DATA/CADDindels.vcf.gz
fi

bedtools intersect -a <(sed '1d' pli.bed) -b patho.vcf.gz | cut -f 5 > tmp/plipatho
bedtools intersect -a <(sed '1d' pli.bed) -b benign.vcf.gz | cut -f 5 > tmp/plibenign

python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p patho.vcf.gz -b benign.vcf.gz 2>/dev/null #2>/dev/null is because the phred score is in the filter column

#paste tmp/ccrbenign tmp/ccrpatho | python ryan/hist.py -o rh.png -t "Clinvar Variant Comparison"
python roccurve.py -t "$TITLE" -c tmp/ccrpatho tmp/ccrbenign -p tmp/plipatho tmp/plibenign -d tmp/caddpatho tmp/caddbenign -g tmp/ccr2patho tmp/ccr2benign
