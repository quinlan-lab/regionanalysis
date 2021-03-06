HOME=/uufs/chpc.utah.edu/common/home/u1021864
DATA=/scratch/ucgd/lustre/u1021864/serial

# making Figure 1

bash runproplot.sh # contains proplot.py

# making Figures 2/5

bash runpathoscore.sh # contains ~/software/pathoscore/pathoscore.py (modified to spit out a pickle dump), stepplot.py, oddsratio.py, combine.py
Rscript genes_with_ccrs_above_pct.R essentials/gnomadbased-ccrs.bed.gz $HOME/public_html/randomplots/genes_with_ccrs_above_pct.pdf

# making Figure 3

bash pfam/pfamvccr.sh # contains pfam/pfamhist.py, pfam/gerppfam.py, pfam/fameval.py
#also outputs supp doc 1 (all pfam histograms) and supp table 3(?) (table contains all stats for each pfam)

# making Figure 4
(printf "chr\tstart\tend\tgene\tccr\tmpcreg\n"; bedtools intersect -a $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz -b essentials/mpc.regions.clean.sorted.bed.gz -wao -sorted | cut -f 1,2,3,4,5,17) | bgzip > ccr.v.mpcregions.bed.gz
Rscript ccr_versus_missensedepletion.R ccr.v.mpcregions.bed.gz 95 $HOME/public_html/randomplots
Rscript ccr_versus_missensedepletion.R ccr.v.mpcregions.bed.gz 99 $HOME/public_html/randomplots
bash constraintcorrelation.sh
# comparison with pLI and missense depletion (NOW Fig 4)
bash venndiagram.sh

# missense variant calculation

#CT=$(python countfuncvars.py essentials/gnomad.vcf.gz) # gives total missense/LoF variants and total (all) variants
#CT=$(echo $CT | cut -f 1 -d " ")
#6748674
cd exacresiduals
bash regions.sh -c -w -v unfiltered -d 10 -d 0
cd -
CT=$(grep VARTRUE exacresiduals/results/unfiltered/exac-regions-novariant.txt | wc -l)
EX=$(awk '{t+=$3-$2} END {print t}' exacresiduals/flatexome.bed)
bc <<< "scale=4; $CT/$EX"

# genes with 99th percentile CCR and how many (Supp Table 1)
Rscript genesw99CCR.R essentials/gnomadbased-ccrs.bed.gz genesw99CCR\(supp_table_1\).tsv

# median lengths and functional variant distance plot
python median.py essentials/gnomadbased-ccrs.bed.gz exacresiduals/results/unfiltered/exac-regions-novariant.txt $HOME/public_html/randomplots/distances.pdf

# gene size vs percentile (Supp Fig 1)
python genevccr.py exacresiduals/flatexome.bed essentials/gnomadbased-ccrs.bed.gz $HOME/public_html/randomplots/genevccr.pdf

# singletons and non-singleton model to estimate saturation, and FDR (Figure 6)

cd exacresiduals
python exac-regions.py -n -w -c "data/exacv2.chr{chrom}.cov.txt.gz" -e data/Homo_sapiens.GRCh37.75.gtf.gz -x data/gnomad-vep-vt.vcf.gz -d 10 -l 0.5 -s data/self-chains.gt90.bed.gz data/segmental.bed.gz -f data/hg19.fa > results/nosingletons10x.5/exac-regions-nosingletons-novariant.txt
python singletons.py -c "data/exacv2.chr{chrom}.cov.txt.gz" -e data/Homo_sapiens.GRCh37.75.gtf.gz -x data/gnomad-vep-vt.vcf.gz -d 10 -l 0.5 -s data/self-chains.gt90.bed.gz data/segmental.bed.gz -f data/hg19.fa | grep -Pv '^X|^Y' > gnomadsingletons.vcf
cd -

SI=$(grep -v "^#" exacresiduals/gnomadsingletons.vcf | wc -l)
TOT=$(zgrep "VARTRUE" essentials/gnomadbased-ccrs.bed.gz | wc -l)
echo "% of singletons in gnomAD\n"
bc <<< "scale=4; $SI/$TOT" 

bash python fdr.sh # contains fdr.py

# code for determining how many genes have no clinvar variation whatsoever but have high CCRs
# used ClinVar variants designated as ANYTHING pathogenic or likely pathogenic
# can't use combine as is, need to extract functional variants and exclude exac

bash clinvarunfiltered.sh
# 95th pct
zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5>=95' | cut -f 4 | sort | uniq > ccr95genes
grep -f clingenes -v -w ccr95genes > notclingenes95
CP=$(wc -l notclingenes95 | cut -d " " -f -1)
CT=$(wc -l ccr95genes | cut -d " " -f -1)
echo "Fraction of genes with CCR >= 95 w/ no known function in ClinVar"
bc <<< "scale=4; ($CP)/$CT" #subtraction means no -v necessary, and because of ccrs in same gene w/o intersection, this works better
zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5>=95' | cut -f 4 | sort | uniq -c | sort -k1,1nr | grep -f clingenes -v -w > notclingenes95ranked
# 99th pct
zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5>=99' | cut -f 4 | sort | uniq > ccr99genes
grep -f clingenes -v -w ccr99genes > notclingenes99
CP=$(wc -l notclingenes99 | cut -d " " -f -1)
CT=$(wc -l ccr99genes | cut -d " " -f -1)
echo "Fraction of genes with CCR >= 99 w/ no known function in ClinVar"
bc <<< "scale=4; ($CP)/$CT" #subtraction means no -v necessary, and because of ccrs in same gene w/o intersection, this works better
zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5>=99' | cut -f 4 | sort | uniq -c | sort -k1,1nr | grep -f clingenes -v -w > notclingenes99ranked

# pfams with no clinvar vars at 99% CCR; list of genes and domains may correlate with EM domains from Kasper's paper
python vars.py -w pathogenic.combine.vcf.gz -f > funcpathos.vcf
zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=99' | bedtools intersect -a stdin -b funcpathos.vcf -v | bedtools intersect -a pfam/pfam.genome.gene.bed.gz -b stdin -sorted -u | bedtools intersect -a stdin -b funcpathos.vcf -v | cut -f 4,5 | sort | uniq -c | sort -k1,1nr | sed 's/^\s*//g' | tr -s " " "\t" > pfamenriched\(supp_table_3\).tsv

# fetal variant comparison
bash fetalvars.sh
bedtools intersect -a <(zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=95') -b fetalvariantsalamillo.vcf fetalvariantscarss.vcf fetalvariantsdrury.vcf

# regression plot for model
python regression.py exacresiduals/results/gnomAD10x.5/weightedresiduals-cpg-novariant.txt $HOME/public_html/randomplots/regression.pdf

# creating C-T synonymous density and all synonymous density v CpG density plots for CCRs >=20 bp, without CpG = 0.
python syndensity.py

# creating C-T synonymous density and all synonymous density v CpG density plots for exons >=20 bp
python cpgexoncalc.py

# pLI v GERP comparison
python genegerp.py

# coverage justification
bash cutofftest.sh

# justification for 20 bp
bash twentycomp.sh

# how many benigns/pathogenics above 95th pct
bash advars.sh
