set -o xtrace

if [ ! -s denovos/mpccontroldenovos.vcf.gz ] | [ ! -s denovos/mpcneurodevdenovos.vcf.gz ]; then
#TODO: add code to make it a proper vcf, vep annotated and everything through varmake.sh
    tr -s "" "\n" < ogfiles/mpcdenovos.txt | sed '1d' | cut -f -5,17 | awk '{print $1,$2,".",$3,$4,"100",$5,"MPC=" $6}' OFS='\t' | sort -k1,1 -k2,2n | grep 'control' | cat mpctrueheader - > denovos/mpccontroldenovos.vcf
    tr -s "" "\n" < ogfiles/mpcdenovos.txt | sed '1d' | cut -f -5,17 | awk '{print $1,$2,".",$3,$4,"100",$5,"MPC=" $6}' OFS='\t' | sort -k1,1 -k2,2n | grep 'id_ddd' | cat mpctrueheader - > denovos/mpcneurodevdenovos.vcf
    bgzip -c denovos/mpccontroldenovos.vcf > denovos/mpccontroldenovos.vcf.gz; tabix denovos/mpccontroldenovos.vcf.gz
    bgzip -c denovos/mpcneurodevdenovos.vcf > denovos/mpcneurodevdenovos.vcf.gz; tabix denovos/mpcneurodevdenovos.vcf.gz
fi

if [ ! -s tmp/neurodev-gnomad.txt ] | [ ! -s tmp/control-gnomad.txt ] ; then
    #bedtools intersect -a denovos/mpcneurodevdenovos.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -wao -sorted > tmp/neurodev-gnomad.txt
    #bedtools intersect -a denovos/mpccontroldenovos.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -wao -sorted > tmp/control-gnomad.txt
    bedtools intersect -a denovos/mpcneurodevdenovos.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -v -header > tmp/neurodev-patho-gnomad.vcf
    bedtools intersect -a denovos/mpccontroldevdenovos.vcf.gz -b $DATA/gnomad-vep-vt.vcf.gz -v -header > tmp/control-benign-gnomad.vcf
    cat <(grep '^#' tmp/neurodev-patho-gnomad.vcf) <(grep -v '^#' tmp/neurodev-patho-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > tmp/neurodev-patho-gnomad.vcf.gz; tabix tmp/neurodev-patho-gnomad.vcf.gz
    cat <(grep '^#' tmp/control-benign-gnomad.vcf) <(grep -v '^#' tmp/control-benign-gnomad.vcf | sort -k1,1 -k2,2n) | bgzip -c > tmp/control-benign-gnomad.vcf.gz; tabix tmp/control-benign-gnomad.vcf.gz
fi

python parvarfilter.py -x tmp/neurodev-gnomad.txt -n neurodev -c -s patho -e gnomad -d genescreens/ad_genecards_clean.txt -f
python parvarfilter.py -x tmp/control-gnomad.txt -n control -c -s benign -e gnomad -d genescreens/ad_genecards_clean.txt -f
cat mpcheader <(sort -k1,1 -k2,2n $DATA/clinvar-patho-gnomad.txt | uniq) > $DATA/clinvar-patho-gnomad.vcf

python scorevars.py -x tmp/neurodev-patho-gnomad.vcf.gz -c gnomad-ccrs.bed.gz -a > tmp/ccr2patho
python scorevars.py -x tmp/control-benign-gnomad.vcf.gz -c gnomad-ccrs.bed.gz -a > tmp/ccr2benign
python caddintersect.py -c $DATA/MPC.vcf.gz -p tmp/neurodev-patho-gnomad.vcf.gz tmp/MPC2patho -b tmp/control-benign-gnomad.vcf.gz tmp/MPC2benign 2>/dev/null #2>/dev/null is because the code is there for CADD and formatting it like a proper VCF isn't important.

paste tmp/MPC2benign tmp/MPC2patho | python hist.py -o mpc_denovo_dist.pdf -t "MPC (based on ExAC v1)" #ExAC v1 based MPC
paste tmp/ccr2benign tmp/ccr2patho | python hist.py -o ccr_denovo_dist.pdf -t "gnomad based CCR" # gnomad based CCR
