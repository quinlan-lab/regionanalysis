# making Figure 1

bash runproplot.sh # contains proplot.py

# making Figures 2/5

bash runpathoscore.sh # contains ~/software/pathoscore/pathoscore.py (modified to spit out a pickle dump), fig2plot.py, oddsratio.py, combine.py

# making Figure 4

bash pfam/pfamvccr.sh # contains pfam/pfamhist.py, pfam/gerppfam.py, pfam/fameval.py

# missense variant calculation

CT=$(python countfuncvars.py essentials/gnomad.vcf.gz) # gives total missense/LoF variants and total (all) variants
CT=$(echo $CT | cut -f 1 -d " ")
#6748674
EX=$(awk '{t+=$3-$2} END {print t}' exacresiduals/flatexome.bed)
bc <<< "scale=4; $CT/$EX"
