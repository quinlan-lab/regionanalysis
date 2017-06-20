#for generating ladder plots and ROC curves, and clinvar filtered files

#--------------------------------------------------------------------------------------
# This is for running the variant annotation script that annotates ClinVar and gnomAD 
# with VEP and decomposes them with VT.  There is also annotation with ExAC and gnomAD
# where relevant.
#--------------------------------------------------------------------------------------

#first runs bash varmake.sh if the clinvar file below is not yet generated

if [ ! -s $DATA/clinvar_20170104-vep-vt.vcf.gz ]; then
    bash varmake.sh $DATA/clinvar_20170104.vcf.gz
fi

if [ ! -s $DATA/gnomad-vep-vt.vcf.gz ]; then
    bash varmake.sh $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz exac # doesn't annotate gnomad with gnomad
fi

#--------------------------------------------------------------------------------------
# This is where we generate the intersections between clinvar, gnomAD, and ExAC for
# filtering out variants down the pipeline.
#--------------------------------------------------------------------------------------

if [ ! -s $DATA/clinvar-exac.txt ] | [ ! -s $DATA/clinvar-gnomad.txt ] ; then
    bedtools intersect -a $DATA/clinvar_20170104-vep-vt.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -wao -sorted > $DATA/clinvar-gnomad.txt
    bedtools intersect -a $DATA/clinvar_20170104-vep-vt.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -wao -sorted > $DATA/clinvar-exac.txt
fi

if [ ! -s $DATA/gnomad-exac.txt ]; then
    bedtools intersect -a $DATA/gnomad-vep-vt.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -wao -sorted > $DATA/gnomad-exac.txt
fi

#--------------------------------------------------------------------------------------
# This is where we generate the metric files necessary for intersecting variants and
# evaluating pLI, RVIS, CADD and CCR.
#--------------------------------------------------------------------------------------

if [ ! -s exac-ccrs.bed.gz ] | [ ! -s gnomad-ccrs.bed.gz ]; then
    sed '1d' exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c > exac-ccrs.bed.gz; tabix exac-ccrs.bed.gz
    sed '1d' exacresiduals/results/newweight30x.5/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c > gnomad-ccrs.bed.gz; tabix gnomad-ccrs.bed.gz
fi

if [ ! -s pli.bed ]; then
    #bedtools groupby -i $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt -g 3,5,6,2 -c 20 -o collapse | tr -s ' ' '\t' > pli.bed
    sed 's/\"//g' exacresiduals/flatexome.bed | sed 's/;//g' > tmp/Homo_sapiens37.bed
    cut -f 2,20 $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt | sed '1d' > tmp/pLI_exac.txt
    awk 'FNR==NR{genes[$1]=$2; next} {for (gene in genes) if (gene == $4) print $0, genes[gene]}' FS='\t' OFS='\t' tmp/pLI_exac.txt tmp/Homo_sapiens37.bed > pli.bed
fi

if [ ! -s rvis.bed ]; then
    sed 's/\"//g' exacresiduals/flatexome.bed | sed 's/;//g' > tmp/Homo_sapiens37.bed
    cut -f 5,11 $DATA/RVIS_Unpublished_ExAC_May2015.txt | sed '1d' > tmp/RVIS_exac.txt
    awk 'FNR==NR{genes[$1]=$2; next} {for (gene in genes) if (gene == $4) print $0, genes[gene]}' FS='\t' OFS='\t' tmp/RVIS_exac.txt tmp/Homo_sapiens37.bed > rvis.bed
fi

if [ ! -s $DATA/CADD.vcf.gz ]; then
    cat caddheader <(sed '1,2d' $DATA/whole_genome_SNVs.tsv | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS="\t") > $DATA/CADD.vcf # whole_genome_SNVs.tsv is the original CADD file
    bgzip -@ 12 -c $DATA/CADD.vcf > $DATA/CADD.vcf.gz; tabix $DATA/CADD.vcf.gz
fi

if [ ! -s $DATA/CADDindels.vcf.gz ]; then
    cat caddheader <(sed '1,2d' $DATA/InDels.tsv | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS="\t") > $DATA/CADDindels.vcf
    bgzip -@ 12 -c $DATA/CADDindels.vcf > $DATA/CADDindels.vcf.gz; tabix $DATA/CADDindels.vcf.gz
fi

if [ ! -s $DATA/fordist_constraint_official_mpc_values.txt.gz ]; then
    #wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/fordist_constraint_official_mpc_values.txt.gz
    #wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/fordist_constraint_official_mpc_values.txt.gz.tbi
    cat mpcheader <(zcat $DATA/fordist_constraint_official_mpc_values.txt.gz | sed '1d' | awk '{print $1,$2,".",$3,$4,"PASS",$NF}' OFS="\t") > $DATA/MPC.vcf #NF because field number is variable
    bgzip $DATA/MPC.vcf -c > $DATA/MPC.vcf.gz; tabix -f $DATA/MPC.vcf.gz
fi

#--------------------------------------------------------------------------------------
# This is where the (CADD + CCR) files are generated.
#--------------------------------------------------------------------------------------

if [ ! -s $DATA/ecCADD.vcf.gz ] | [ ! -s $DATA/ecCADDindels.vcf.gz ]; then
    cat <(zgrep "^#" $DATA/CADD.vcf.gz) <(bedtools intersect -a $DATA/CADD.vcf.gz -b exac-ccrs.bed.gz -wb -sorted | awk '{print $1,$2,$3,$4,$5,$6,$7+-10*log(1-$21/100)/log(10)}' FS='\t' OFS='\t') > $DATA/ecCADD.vcf
    bgzip -@ 12 -c $DATA/ecCADD.vcf > $DATA/ecCADD.vcf.gz; tabix $DATA/ecCADD.vcf.gz
    cat <(zgrep "^#" $DATA/CADDindels.vcf.gz) <(bedtools intersect -a $DATA/CADDindels.vcf.gz -b exac-ccrs.bed.gz -wb -sorted | awk '{print $1,$2,$3,$4,$5,$6,$7+-10*log(1-$21/100)/log(10)}' FS='\t' OFS='\t') > $DATA/ecCADDindels.vcf
    bgzip -@ 12 -c $DATA/ecCADDindels.vcf > $DATA/ecCADDindels.vcf.gz; tabix $DATA/ecCADDindels.vcf.gz
fi

if [ ! -s $DATA/gcCADD.vcf.gz ] | [ ! -s $DATA/gcCADDindels.vcf.gz ]; then
    cat <(zgrep "^#" $DATA/CADD.vcf.gz) <(bedtools intersect -a $DATA/CADD.vcf.gz -b gnomad-ccrs.bed.gz -wb -sorted | awk '{print $1,$2,$3,$4,$5,$6,$7+-10*log(1-$21/100)/log(10)}' FS='\t' OFS='\t') > $DATA/gcCADD.vcf
    bgzip -@ 12 -c $DATA/gcCADD.vcf > $DATA/gcCADD.vcf.gz; tabix $DATA/gcCADD.vcf.gz
    cat <(zgrep "^#" $DATA/CADDindels.vcf.gz) <(bedtools intersect -a $DATA/CADDindels.vcf.gz -b gnomad-ccrs.bed.gz -wb -sorted | awk '{print $1,$2,$3,$4,$5,$6,$7+-10*log(1-$21/100)/log(10)}' FS='\t' OFS='\t') > $DATA/gcCADDindels.vcf
    bgzip -@ 12 -c $DATA/gcCADDindels.vcf > $DATA/gcCADDindels.vcf.gz; tabix $DATA/gcCADDindels.vcf.gz
fi

#--------------------------------------------------------------------------------------
# Here, we generate the variant files filtered on ExAC and gnomAD and also store how
# many variants there are in total.
#--------------------------------------------------------------------------------------

# generates the "patho.vcf" and "benign.vcf" files that are strictly filtered based on our criteria
python parvarfilter.py -x $DATA/clinvar-gnomad.txt -n clinvar -c -s patho -e gnomad -d genescreens/ad_genecards_clean.txt -f #-i genescreens/clingen_level3_genes_2015_02_27.tsv
cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-patho-gnomad.txt | uniq) > $DATA/clinvar-patho-gnomad.vcf
python parvarfilter.py -x $DATA/clinvar-exac.txt -n clinvar -c -s patho -e exac -d genescreens/ad_genecards_clean.txt -f #-i genescreens/clingen_level3_genes_2015_02_27.tsv
cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-patho-exac.txt | uniq) > $DATA/clinvar-patho-exac.vcf
python parvarfilter.py -x $DATA/clinvar-gnomad.txt -n clinvar -c -s benign -e gnomad -d genescreens/ad_genecards_clean.txt -f #-i genescreens/clingen_level3_genes_2015_02_27.tsv
cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-benign-gnomad.txt | uniq) > $DATA/clinvar-benign-gnomad.vcf

EP=$(grep -v "^#" $DATA/clinvar-patho-exac.vcf | wc -l)
GP=$(grep -v "^#" $DATA/clinvar-patho-gnomad.vcf | wc -l)
GB=$(grep -v "^#" $DATA/clinvar-benign-gnomad.vcf | wc -l)

#--------------------------------------------------------------------------------------
# This is where we create the score files for CCR based on a weighted average (in case
# of INDELs) and also generate which variants successfully overlap with CCR to create
# a common set of variants across all metrics.  Except for ExAC-based CCR benign vars,
# which are in the getopts script below.
#--------------------------------------------------------------------------------------

#exac
cat <(grep '^#' $DATA/clinvar-patho-exac.vcf) <(grep -v '^#' $DATA/clinvar-patho-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-patho-exac.vcf.gz; tabix $DATA/clinvar-patho-exac.vcf.gz
python scorevars.py -x $DATA/clinvar-patho-exac.vcf.gz -c exac-ccrs.bed.gz -a > tmp/ccrpatho

#gnomAD
cat <(grep '^#' $DATA/clinvar-benign-gnomad.vcf) <(grep -v '^#' $DATA/clinvar-benign-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-benign-gnomad.vcf.gz; tabix $DATA/clinvar-benign-gnomad.vcf.gz
cat <(grep '^#' $DATA/clinvar-patho-gnomad.vcf) <(grep -v '^#' $DATA/clinvar-patho-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-patho-gnomad.vcf.gz; tabix $DATA/clinvar-patho-gnomad.vcf.gz
python scorevars.py -x $DATA/clinvar-patho-gnomad.vcf.gz -c exac-ccrs.bed.gz -a > tmp/ccr2patho
python scorevars.py -x $DATA/clinvar-benign-gnomad.vcf.gz -c exac-ccrs.bed.gz -a > tmp/ccr2benign

#--------------------------------------------------------------------------------------
# This is where we create the score files for benign variants (depending on ClinVar
# choice or gnomAD benign choice).  Additionally CADD pathogenic scores are generated
# here because it is way faster to work with CADD two files at a time than doing two
# passes.
#--------------------------------------------------------------------------------------

while getopts ":t:g:c" opt; do
    case $opt in
        t)
            TITLE+=("$OPTARG")
            ;;
        g)
            AF=$OPTARG
            echo "-gnomad as benign input was triggered" >&2
            if [ ! -s $DATA/gnomad-benign-exac.vcf ]; then
                python parvarfilter.py -x $DATA/gnomad-exac.txt -n gnomad -s benign -e exac -d genescreens/ad_genecards_clean.txt -f #-i genescreens/clingen_level3_genes_2015_02_27.tsv
                cat <(zgrep "^#" $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/gnomad-benign-exac.txt | uniq) > $DATA/gnomad-benign-exac.vcf
            fi
            python secondfilter.py -x $DATA/gnomad-benign-exac.vcf -f $AF -d genescreens/ad_genecards_clean.txt

            #exac
            cat <(grep '^#' $DATA/gnomad-benign-exac-filtered.vcf) <(grep -v '^#' $DATA/gnomad-benign-exac-filtered.vcf | shuf -n $EP | sort -k1,1 -k2,2n) | bgzip -c > $DATA/gnomad-benign-exac.vcf.gz; tabix $DATA/gnomad-benign-exac.vcf.gz # shuf adds subsampling
            EB=$(zgrep -v "^#" $DATA/gnomad-benign-exac.vcf.gz | wc -l)
            #bedtools intersect -a $DATA/gnomad-benign-exac.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -v > $DATA/gnomadbenigns.vcf.gz
            python scorevars.py -x $DATA/gnomad-benign-exac.vcf.gz -c exac-ccrs.bed.gz -a > tmp/ccrbenign
            python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/gnomad-benign-exac.vcf.gz -f tmp/caddpatho tmp/caddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column
            python caddintersect.py -c $DATA/ecCADD.vcf.gz -d $DATA/ecCADDindels.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/gnomad-benign-exac.vcf.gz -f tmp/eccaddpatho tmp/eccaddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column so cyvcf2 spits an error that is ignorable over and over
            python caddintersect.py -c $DATA/MPC.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/gnomad-benign-exac.vcf.gz -f tmp/MPCpatho tmp/MPCbenign 2>/dev/null #2>/dev/null is because the code is there for CADD and formatting it like a proper VCF isn't important.
            bedtools intersect -a rvis.bed -b $DATA/gnomad-benign-exac.vcf.gz | cut -f 5 > tmp/rvisbenign
            bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/gnomad-benign-exac.vcf.gz | cut -f 5 > tmp/plibenign
            ;;
        c)
            echo "-clinvar input triggered" >&2
            python parvarfilter.py -x $DATA/clinvar-exac.txt -n clinvar -c -s benign -e exac -d genescreens/ad_genecards_clean.txt -f #-i genescreens/clingen_level3_genes_2015_02_27.tsv
            cat <(zgrep "^#" $DATA/clinvar_20170104.vcf.gz ) <(sort -k1,1 -k2,2n $DATA/clinvar-benign-exac.txt | uniq) > $DATA/clinvar-benign-exac.vcf
            EB=$(grep -v "^#" $DATA/clinvar-benign-exac.vcf | wc -l)

            #exac
            cat <(grep '^#' $DATA/clinvar-benign-exac.vcf) <(grep -v '^#' $DATA/clinvar-benign-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/clinvar-benign-exac.vcf.gz; tabix $DATA/clinvar-benign-exac.vcf.gz
            python scorevars.py -x $DATA/clinvar-benign-exac.vcf.gz -c exac-ccrs.bed.gz -a > tmp/ccrbenign
            python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/clinvar-benign-exac.vcf.gz -f tmp/caddpatho tmp/caddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column
            python caddintersect.py -c $DATA/ecCADD.vcf.gz -d $DATA/ecCADDindels.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/clinvar-benign-exac.vcf.gz -f tmp/eccaddpatho tmp/eccaddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column so cyvcf2 spits an error that is ignorable over and over
            python caddintersect.py -c $DATA/MPC.vcf.gz -p $DATA/clinvar-patho-exac.vcf.gz -b $DATA/clinvar-benign-exac.vcf.gz -f tmp/MPCpatho tmp/MPCbenign 2>/dev/null #2>/dev/null is because the code is there for CADD and formatting it like a proper VCF isn't important.
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

#--------------------------------------------------------------------------------------
# This is where the variant files for McRae (DDD study) and deLigt (ID study) are
# generated (deLigt is WIP atm)
#--------------------------------------------------------------------------------------

echo "-title input: '${TITLE[@]}'" >&2

if [ ! -s $DATA/mcrae-patho-exac.vcf.gz ] | [ ! -s $DATA/mcrae-patho-gnomad.vcf.gz ]; then
    cat vcfheader <(sed '1,2d' denovos/mcrae.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/mcrae.vcf.gz; tabix denovos/mcrae.vcf.gz
    bash varmake.sh denovos/mcrae.vcf.gz
    bedtools intersect -a $DATA/mcrae-vep-vt.vcf.gz -b $DATA/ExAC.r1.vt.vep.vcf.gz -wao > $DATA/mcrae-exac.txt #-sorted
    bedtools intersect -a $DATA/mcrae-vep-vt.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -wao > $DATA/mcrae-gnomad.txt
    python parvarfilter.py -x $DATA/mcrae-exac.txt -e exac -f -n mcrae -s patho
    python parvarfilter.py -x $DATA/mcrae-gnomad.txt -e gnomad -f -n mcrae -s patho
    cat <(grep '^#' $DATA/mcrae-patho-exac.vcf) <(grep -v '^#' $DATA/mcrae-patho-exac.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/mcrae-patho-exac.vcf.gz; tabix $DATA/mcrae-patho-exac.vcf.gz
    cat <(grep '^#' $DATA/mcrae-patho-gnomad.vcf) <(grep -v '^#' $DATA/mcrae-patho-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > $DATA/mcrae-patho-gnomad.vcf.gz; tabix $DATA/mcrae-patho-gnomad.vcf.gz
fi

if [ ! -s denovos/deligt.vcf.gz ] || [ ! -s denovos/deligtcontrol.vcf.gz ]; then
    cat vcfheader <(sed '1,2d' denovos/deligt.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/deligt.vcf.gz; tabix denovos/deligt.vcf.gz
    cat vcfheader <(sed '1,2d' denovos/deligtcontrol.bed | awk '{print $1,$3,".",$4,$5,".","PASS","GENEINFO=" $6 ";"}' OFS="\t" | sort -k1,1 -k2,2n) | bgzip -c > denovos/deligtcontrol.vcf.gz; tabix denovos/deligtcontrol.vcf.gz
fi

#--------------------------------------------------------------------------------------
# This is where the score files for McRae (DDD study) and deLigt (ID study) are
# generated (deLigt is WIP atm)
#--------------------------------------------------------------------------------------

bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-patho-exac.vcf.gz | cut -f 5 > tmp/plipatho
bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-patho-gnomad.vcf.gz | cut -f 5 > tmp/pli2patho
bedtools intersect -a <(sed '1d' pli.bed) -b $DATA/clinvar-benign-gnomad.vcf.gz | cut -f 5 > tmp/pli2benign

bedtools intersect -a rvis.bed -b $DATA/clinvar-patho-exac.vcf.gz | cut -f 5 > tmp/rvispatho
bedtools intersect -a rvis.bed -b $DATA/clinvar-patho-gnomad.vcf.gz | cut -f 5 > tmp/rvis2patho
bedtools intersect -a rvis.bed -b $DATA/clinvar-benign-gnomad.vcf.gz | cut -f 5 > tmp/rvis2benign

python scorevars.py -x $DATA/mcrae-patho-exac.vcf.gz -c exac-ccrs.bed.gz -a > tmp/mcraepatho
python scorevars.py -x $DATA/mcrae-patho-gnomad.vcf.gz -c exac-ccrs.bed.gz -a > tmp/mcrae2patho

python caddintersect.py -c $DATA/CADD.vcf.gz -d $DATA/CADDindels.vcf.gz -p $DATA/clinvar-patho-gnomad.vcf.gz -b $DATA/clinvar-benign-gnomad.vcf.gz -f tmp/cadd2patho tmp/cadd2benign 2>/dev/null #2>/dev/null is because the phred score is in the filter column so cyvcf2 spits an error that is ignorable over and over

python caddintersect.py -c $DATA/gcCADD.vcf.gz -d $DATA/gcCADDindels.vcf.gz -p $DATA/clinvar-patho-gnomad.vcf.gz -b $DATA/clinvar-benign-gnomad.vcf.gz -f tmp/gccaddpatho tmp/gccaddbenign 2>/dev/null #2>/dev/null is because the phred score is in the filter column so cyvcf2 spits an error that is ignorable over and over

python caddintersect.py -c $DATA/MPC.vcf.gz -p $DATA/clinvar-patho-gnomad.vcf.gz -b $DATA/clinvar-benign-gnomad.vcf.gz -f tmp/MPC2patho tmp/MPC2benign 2>/dev/null #2>/dev/null is because the code is there for CADD and formatting it like a proper VCF isn't important.

paste tmp/mcraepatho | python hist.py -o mcraeexac_dist.pdf #-t "$TITLE"
paste tmp/mcrae2patho | python hist.py -o mcraegnomad_dist.pdf #-t "$TITLE"
paste tmp/ccrbenign tmp/ccrpatho | python hist.py -o exac_dist.pdf -t "${TITLE[0]}"
paste tmp/ccr2benign tmp/ccr2patho | python hist.py -o gnomad_dist.pdf -t "${TITLE[1]}"

python roccurve.py -t "${TITLE[0]}" -c tmp/ccrpatho tmp/ccrbenign -p tmp/plipatho tmp/plibenign -d tmp/caddpatho tmp/caddbenign -r tmp/rvispatho tmp/rvisbenign -b tmp/eccaddpatho tmp/eccaddbenign -m tmp/MPCpatho tmp/MPCbenign -o exac_roc.pdf -n $EP $EB
python roccurve.py -t "${TITLE[1]}" -g tmp/ccr2patho tmp/ccr2benign -p tmp/pli2patho tmp/pli2benign -d tmp/cadd2patho tmp/cadd2benign -r tmp/rvis2patho tmp/rvis2benign -b tmp/gccaddpatho tmp/gccaddbenign -m tmp/MPC2patho tmp/MPC2benign -o gnomad_roc.pdf -n $GP $GB
