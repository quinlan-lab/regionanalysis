# making Figure 1

bash runproplot.sh # contains proplot.py

# making Figures 2/5

bash runpathoscore.sh # contains ~/software/pathoscore/pathoscore.py (modified to spit out a pickle dump), fig2plot.py, oddsratio.py, combine.py
Rscript genes_with_ccrs_above_pct.R

# making Figure 3 (aaron has code for Figure 4)

bash pfam/pfamvccr.sh # contains pfam/pfamhist.py, pfam/gerppfam.py, pfam/fameval.py

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

# median lengths and functional variant distance plot
python median.py /uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/distances.pdf

# singletons and non-singleton model to estimate FDR (Figure 6)

cd exacresiduals
python exac-regions.py -n -w -c "data/exacv2.chr{chrom}.cov.txt.gz" -e data/Homo_sapiens.GRCh37.75.gtf.gz -x data/gnomad-vep-vt.vcf.gz -d 10 -l 0.5 -s data/self-chains.gt90.bed.gz data/segmental.bed.gz -f data/hg19.fa > results/nosingletons10x.5/exac-regions-nosingletons-novariant.txt
python singletons.py -c "data/exacv2.chr{chrom}.cov.txt.gz" -e data/Homo_sapiens.GRCh37.75.gtf.gz -x data/gnomad-vep-vt.vcf.gz -d 10 -l 0.5 -s data/self-chains.gt90.bed.gz data/segmental.bed.gz -f data/hg19.fa | grep -Pv '^X|^Y' > gnomadsingletons.vcf
cd -

SI=$(grep -v "^#" exacresiduals/gnomadsingletons.vcf | wc -l)
TOT=$(zgrep "VARTRUE" essentials/gnomadbased-ccrs.bed.gz | wc -l)
echo "% of singletons in gnomAD\n"
bc <<< "scale=4; $SI/$TOT" 

bash python fdr.sh # contains fdr.py

# code for determining how many genes have no clinvar variation but have high CCRs
# used ClinVar variants designated as functional and not overlapping ExAC
# can't use combine as is, need to extract functional variants and exclude exac

python vars.py -w pathogenic.combine.vcf.gz -f > funcpathos.vcf
CP=$(zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=95' | bedtools intersect -a stdin -b funcpathos.vcf | cut -f 4 | sort | uniq | wc -l)
CT=$(zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=95' | cut -f 4 | sort | uniq | wc -l)
echo "Fraction of genes with CCR >= 95 w/ no known function in ClinVar\n"
bc <<< "scale=4; ($CT-$CP)/$CT"
CP=$(zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=99' | bedtools intersect -a stdin -b funcpathos.vcf | cut -f 4 | sort | uniq | wc -l)
CT=$(zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=99' | cut -f 4 | sort | uniq | wc -l)
echo "Fraction of genes with CCR >= 99 w/ no known function in ClinVar\n"
bc <<< "scale=4; ($CT-$CP)/$CT"
