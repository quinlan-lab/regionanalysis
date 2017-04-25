#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=6:00:00
#set -exo pipefail -o nounset

#for generating ladder plots and ROC curves
#first run bash clinvarmake.sh if the clinvar file below is not yet generated
TITLE=$1

python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e exac -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s patho #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e gnomad -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s patho #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e exac -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s benign #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e gnomad -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s benign #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria

#exac
cat <(grep '^#' clinvar-benign-exac.vcf) <(grep -v '^#' clinvar-benign-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > clinvar-benign-exac.vcf.gz; tabix clinvar-benign-exac.vcf.gz
cat <(grep '^#' clinvar-patho-exac.vcf) <(grep -v '^#' clinvar-patho-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > clinvar-patho-exac.vcf.gz; tabix clinvar-patho-exac.vcf.gz
bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b clinvar-patho-exac.vcf.gz | cut -f 14 > tmp/ccrpatho
bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b clinvar-benign-exac.vcf.gz | cut -f 14 > tmp/ccrbenign

#gnomAD
cat <(grep '^#' clinvar-benign-gnomad.vcf) <(grep -v '^#' clinvar-benign-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > clinvar-benign-gnomad.vcf.gz; tabix clinvar-benign-gnomad.vcf.gz
cat <(grep '^#' clinvar-patho-gnomad.vcf) <(grep -v '^#' clinvar-patho-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > clinvar-patho-gnomad.vcf.gz; tabix clinvar-patho-gnomad.vcf.gz
bedtools intersect -a <(sed '1d' exacresiduals/results/newweight30x.5/weightedresiduals-cpg-novariant.txt) -b clinvar-patho-gnomad.vcf.gz | cut -f 14 > tmp/ccr2patho
bedtools intersect -a <(sed '1d' exacresiduals/results/newweight30x.5/weightedresiduals-cpg-novariant.txt) -b clinvar-benign-gnomad.vcf.gz | cut -f 14 > tmp/ccr2benign

# generate bed file from pli file (and a vcf CADD file by manually adding a header to the 1.3 version TSV file)

if [ ! -s pli.bed ]; then
    bedtools groupby -i $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt -g 3,5,6,2 -c 20 -o collapse | tr -s ' ' '\t' > pli.bed
fi

if [ ! -s rvis.bed ]; then
    sed 's/\"//g' exacresiduals/flatexome.bed | sed 's/;//g' > tmp/Homo_sapiens37.bed
    cut -f 5,11 $DATA/RVIS_Unpublished_ExAC_May2015.txt | sed '1d' > tmp/RVIS_exac.txt
    awk 'FNR==NR{genes[$1]=$2; next} {for (gene in genes) if (gene == $4) print $0, genes[gene]}' FS='\t' OFS='\t' tmp/RVIS_exac.txt tmp/Homo_sapiens37.bed > rvis.bed
fi

if [ ! -s denovos/mcrae.vcf.gz ]; then
    cat vcfheader <(sed '1,2d' denovos/mcrae.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/mcrae.vcf.gz; tabix denovos/mcrae.vcf.gz
    bash varmake.sh denovos/mcrae.vcf.gz
    python varfilter.py -x $DATA/mcrae-vep-anno-vt.vcf.gz -e exac -f -n mcrae
fi

if [ ! -s denovos/deligt.vcf.gz ] || [ ! -s denovos/deligtcontrol.vcf.gz ]; then
    cat vcfheader <(sed '1,2d' denovos/deligt.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/deligt.vcf.gz; tabix denovos/deligt.vcf.gz
    cat vcfheader <(sed '1,2d' denovos/deligtcontrol.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/deligtcontrol.vcf.gz; tabix denovos/deligtcontrol.vcf.gz
fi

if [ ! -s $DATA/CADD.vcf.gz ]; then
    cat caddheader <(sed '1,2d' $DATA/whole_genome_SNVs.tsv | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS="\t") > $DATA/CADD.vcf # whole_genome_SNVs.tsv is the original CADD file
    bgzip -@ 12 -c $DATA/CADD.vcf > $DATA/CADD.vcf.gz; tabix $DATA/CADD.vcf.gz
fi

if [ ! -s $DATA/CADDindels.vcf.gz ]; then
    cat caddheader <(sed '1,2d' $DATA/InDels.tsv | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS="\t") > $DATA/CADDindels.vcf
    bgzip -@ 12 -c $DATA/CADDindels.vcf > $DATA/CADDindels.vcf.gz; tabix $DATA/CADDindels.vcf.gz
fi

bedtools intersect -a <(sed '1d' pli.bed) -b patho-exac.vcf.gz | cut -f 5 > tmp/plipatho
bedtools intersect -a <(sed '1d' pli.bed) -b benign-exac.vcf.gz | cut -f 5 > tmp/plibenign

bedtools intersect -a rvis.bed -b patho-exac.vcf.gz | cut -f 5 > tmp/rvispatho
bedtools intersect -a rvis.bed -b benign-exac.vcf.gz | cut -f 5 > tmp/rvisbenign

bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b mcrae-exac.vcf | cut -f 14 > tmp/mcraepatho

python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p patho-gnomad.vcf.gz -b benign-gnomad.vcf.gz 2>/dev/null #2>/dev/null is because the phred score is in the filter column # don't forget to switch gnomad for exac and vice-versa

paste tmp/mcraepatho | python hist.py -o mcrae_dist.pdf -t "$TITLE"
paste tmp/ccrbenign tmp/ccrpatho | python hist.py -o exac_dist.pdf -t "$TITLE"
paste tmp/ccr2benign tmp/ccr2patho | python hist.py -o gnomad_dist.pdf -t "$TITLE"
python roccurve.py -t "$TITLE" -c tmp/ccrpatho tmp/ccrbenign -p tmp/plipatho tmp/plibenign -d tmp/caddpatho tmp/caddbenign -g tmp/ccr2patho tmp/ccr2benign -r tmp/rvispatho tmp/rvisbenign
