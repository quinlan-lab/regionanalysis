#for generating ladder plots and ROC curves
#first run bash clinvarmake.sh if the clinvar file below is not yet generated

if [ ! -s $DATA/clinvar_20170104-vep-vt.vcf.gz ]; then
    bash varmake.sh $DATA/clinvar_20170104.vcf.gz
fi

if [ ! -s $DATA/gnomad-vep-vt.vcf.gz ]; then
    bash varmake.sh $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz exac # doesn't annotate gnomad with gnomad
fi

if [ ! -s $DATA/clinvar-exac.txt ] | [ ! -s $DATA/clinvar-gnomad.txt ] ; then
    bedtools intersect -a $DATA/clinvar_20170104-vep-vt.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -wa -wb > $DATA/clinvar-gnomad.txt
    bedtools intersect -a $DATA/clinvar_20170104-vep-vt.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -wa -wb > $DATA/clinvar-exac.txt
fi

if [ ! -s $DATA/gnomad-exac.txt ]; then
    bedtools intersect -a $DATA/gnomad-vep-vt.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -wa -wb > $DATA/gnomad-exac.txt
fi

python parvarfilter.py -x $DATA/clinvar-gnomad.txt -f -n clinvar -c -s patho -e gnomad -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv
cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-patho-gnomad.txt | uniq) > $DATA/clinvar-patho-gnomad.vcf
python parvarfilter.py -x $DATA/clinvar-exac.txt -f -n clinvar -c -s patho -e exac -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv
cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-patho-exac.txt | uniq) > $DATA/clinvar-patho-exac.vcf
python parvarfilter.py -x $DATA/clinvar-gnomad.txt -f -n clinvar -c -s benign -e gnomad -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv
cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-benign-gnomad.txt | uniq) > $DATA/clinvar-benign-gnomad.vcf
#python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e exac -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s patho #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
#python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e gnomad -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s patho #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
#python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e gnomad -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s benign #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria

EP=$(grep -v "^#" $DATA/clinvar-patho-exac.vcf | wc -l)
GP=$(grep -v "^#" $DATA/clinvar-patho-gnomad.vcf | wc -l)
GB=$(grep -v "^#" $DATA/clinvar-benign-gnomad.vcf | wc -l)

#exac
cat <(grep '^#' $DATA/clinvar-patho-exac.vcf) <(grep -v '^#' $DATA/clinvar-patho-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-patho-exac.vcf.gz; tabix $DATA/clinvar-patho-exac.vcf.gz
bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b $DATA/clinvar-patho-exac.vcf.gz | cut -f 14 > tmp/ccrpatho

#gnomAD
cat <(grep '^#' $DATA/clinvar-benign-gnomad.vcf) <(grep -v '^#' $DATA/clinvar-benign-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-benign-gnomad.vcf.gz; tabix $DATA/clinvar-benign-gnomad.vcf.gz
cat <(grep '^#' $DATA/clinvar-patho-gnomad.vcf) <(grep -v '^#' $DATA/clinvar-patho-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-patho-gnomad.vcf.gz; tabix $DATA/clinvar-patho-gnomad.vcf.gz
bedtools intersect -a <(sed '1d' exacresiduals/results/newweight30x.5/weightedresiduals-cpg-novariant.txt) -b $DATA/clinvar-patho-gnomad.vcf.gz | cut -f 14 > tmp/ccr2patho
bedtools intersect -a <(sed '1d' exacresiduals/results/newweight30x.5/weightedresiduals-cpg-novariant.txt) -b $DATA/clinvar-benign-gnomad.vcf.gz | cut -f 14 > tmp/ccr2benign

while getopts ":t:gc" opt; do
    case $opt in
        t)
            TITLE+=("$OPTARG")
            ;;
        g)
            echo "-gnomad as benign input was triggered" >&2
            python parvarfilter.py -x $DATA/gnomad-exac.txt -f -n gnomad -c -s benign -e exac -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv
            cat <(zgrep "^#" $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/gnomad-benign-exac.txt | uniq) > $DATA/gnomad-benign-exac.vcf
            #python varfilter.py -x $DATA/gnomad-vep-anno-vt.vcf.gz -e exac -f -n gnomad -s benign -d genescreens/ad_genecards_clean.txt #ogfiles/all_ad.tsv # -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria

            #exac
            cat <(grep '^#' $DATA/gnomad-benign-exac.vcf) <(grep -v '^#' $DATA/gnomad-benign-exac.vcf | shuf -n $EP | sort -k1,1 -k2,2n) | bgzip -c > $DATA/gnomad-benign-exac.vcf.gz; tabix $DATA/gnomad-benign-exac.vcf.gz # shuf adds subsampling
            EB=$(zgrep -v "^#" $DATA/gnomad-benign-exac.vcf.gz | wc -l)
            #bedtools intersect -a $DATA/gnomad-benign-exac.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -v > $DATA/gnomadbenigns.vcf.gz
            bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b $DATA/gnomad-benign-exac.vcf.gz | cut -f 14 > tmp/ccrbenign
            python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/gnomad-benign-exac.vcf.gz -f tmp/caddpatho tmp/caddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column
            bedtools intersect -a rvis.bed -b $DATA/gnomad-benign-exac.vcf.gz | cut -f 5 > tmp/rvisbenign
            bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/gnomad-benign-exac.vcf.gz | cut -f 5 > tmp/plibenign
            ;;
        c)
            echo "-clinvar input triggered" >&2
            python parvarfilter.py -x $DATA/clinvar-exac.txt -f -n clinvar -c -s benign -e exac -d genescreens/ad_genecards_clean.txt #-i genescreens/clingen_level3_genes_2015_02_27.tsv
            cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-benign-exac.txt | uniq) > $DATA/clinvar-benign-exac.vcf
            #python varfilter.py -x $DATA/clinvar_20170104-vep-anno-vt.vcf.gz -e exac -d genescreens/ad_genecards_clean.txt -c -f -n clinvar -s benign #-i genescreens/clingen_level3_genes_2015_02_27.tsv # generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
            EB=$(grep -v "^#" $DATA/clinvar-benign-exac.vcf | wc -l)

            #exac
            cat <(grep '^#' $DATA/clinvar-benign-exac.vcf) <(grep -v '^#' $DATA/clinvar-benign-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-benign-exac.vcf.gz; tabix $DATA/clinvar-benign-exac.vcf.gz
            bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b $DATA/clinvar-benign-exac.vcf.gz | cut -f 14 > tmp/ccrbenign
            python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/clinvar-benign-exac.vcf.gz -f tmp/caddpatho tmp/caddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column
            bedtools intersect -a rvis.bed -b $DATA/clinvar-benign-exac.vcf.gz | cut -f 5 > tmp/rvisbenign
            bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-benign-exac.vcf.gz | cut -f 5 > tmp/plibenign
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

echo "-title input: '${TITLE[@]}'" >&2

if [ ! -s pli.bed ]; then
    bedtools groupby -i $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt -g 3,5,6,2 -c 20 -o collapse | tr -s ' ' '\t' > pli.bed
fi

if [ ! -s rvis.bed ]; then
    sed 's/\"//g' exacresiduals/flatexome.bed | sed 's/;//g' > tmp/Homo_sapiens37.bed
    cut -f 5,11 $DATA/RVIS_Unpublished_ExAC_May2015.txt | sed '1d' > tmp/RVIS_exac.txt
    awk 'FNR==NR{genes[$1]=$2; next} {for (gene in genes) if (gene == $4) print $0, genes[gene]}' FS='\t' OFS='\t' tmp/RVIS_exac.txt tmp/Homo_sapiens37.bed > rvis.bed
fi

if [ ! -s $DATA/mcrae-vep-anno-vt.vcf.gz ]; then
    cat vcfheader <(sed '1,2d' denovos/mcrae.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/mcrae.vcf.gz; tabix denovos/mcrae.vcf.gz
    bash varmake.sh denovos/mcrae.vcf.gz
    python varfilter.py -x $DATA/mcrae-vep-anno-vt.vcf.gz -e exac -f -n mcrae -s patho
    python varfilter.py -x $DATA/mcrae-vep-anno-vt.vcf.gz -e gnomad -f -n mcrae -s patho
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

bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-patho-exac.vcf.gz | cut -f 5 > tmp/plipatho
bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-patho-gnomad.vcf.gz | cut -f 5 > tmp/pli2patho
bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-benign-gnomad.vcf.gz | cut -f 5 > tmp/pli2benign

bedtools intersect -a rvis.bed -b $DATA/clinvar-patho-exac.vcf.gz | cut -f 5 > tmp/rvispatho
bedtools intersect -a rvis.bed -b $DATA/clinvar-patho-gnomad.vcf.gz | cut -f 5 > tmp/rvis2patho
bedtools intersect -a rvis.bed -b $DATA/clinvar-benign-gnomad.vcf.gz | cut -f 5 > tmp/rvis2benign

bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b $DATA/mcrae-patho-exac.vcf | cut -f 14 > tmp/mcraepatho
bedtools intersect -a <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -b $DATA/mcrae-patho-gnomad.vcf | cut -f 14 > tmp/mcrae2patho

python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p $DATA/clinvar-patho-gnomad.vcf.gz -b $DATA/clinvar-benign-gnomad.vcf.gz -f tmp/cadd2patho tmp/cadd2benign 2>/dev/null #2>/dev/null is because the phred score is in the filter column so cyvcf2 spits an error that is ignorable over and over # don't forget to switch gnomad for exac and vice-versa

paste tmp/mcraepatho | python hist.py -o mcraeexac_dist.pdf #-t "$TITLE"
paste tmp/mcrae2patho | python hist.py -o mcraegnomad_dist.pdf #-t "$TITLE"
paste tmp/ccrbenign tmp/ccrpatho | python hist.py -o exac_dist.pdf -t "${TITLE[0]}"
paste tmp/ccr2benign tmp/ccr2patho | python hist.py -o gnomad_dist.pdf -t "${TITLE[1]}"
python roccurve.py -t "${TITLE[0]}" -c tmp/ccrpatho tmp/ccrbenign -p tmp/plipatho tmp/plibenign -d tmp/caddpatho tmp/caddbenign -r tmp/rvispatho tmp/rvisbenign -o exac_roc.pdf -n $EP $EB
python roccurve.py -t "${TITLE[1]}" -g tmp/ccr2patho tmp/ccr2benign -p tmp/pli2patho tmp/pli2benign -d tmp/cadd2patho tmp/cadd2benign -r tmp/rvis2patho tmp/rvis2benign -o gnomad_roc.pdf -n $GP $GB

if [ ! -s $DATA/newexac.vcf ]; then
    bedtools intersect -a $DATA/CADD.vcf.gz -b <(sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt) -wb | awk '{print $1,$2,$3,$4,$5,$6,$7+-10*log(1-$21/100)/log(10)}' FS='\t' OFS='\t' > $DATA/newexac.vcf
fi
