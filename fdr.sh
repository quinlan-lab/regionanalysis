# estimate saturation by comparing singleton based CCR to doubleton or > based CCR #
sed '1d' exacresiduals/results/gnomAD10x.5/weightedresiduals-cpg-novariant.txt | awk '{n=split($7,a,/,/); split(a[1],s,/-/); m=split(a[n],e,/-/); {printf "%s\t%s\t%s", $1, s[1], e[m]} {for (x=4; x<=NF; x++) printf "\t%s", $x} {printf "\n"}}' OFS="\t" | uniq > tmp/ccrsingenomespace.bed
sed '1d' exacresiduals/results/nosingletons10x.5/weightedresiduals-cpg-nosingletons-novariant.txt | awk '{n=split($7,a,/,/); split(a[1],s,/-/); m=split(a[n],e,/-/); {printf "%s\t%s\t%s", $1, s[1], e[m]} {for (x=4; x<=NF; x++) printf "\t%s", $x} {printf "\n"}}' OFS="\t" | uniq > tmp/doubletonccrsingenomespace.bed
bedtools intersect -a tmp/doubletonccrsingenomespace.bed -b tmp/ccrsingenomespace.bed -wa -wb > ccrintersect.txt
#bedtools intersect -a <(sed '1d' exacresiduals/results/gnomAD10x.5/weightedresiduals-cpg-novariant.txt) -b <(sed '1d' exacresiduals/results/nosingletons10x.5/weightedresiduals-cpg-nosingletons-novariant.txt) -wa -wb > ccrintersect.txt
python saturation.py ccrintersect.txt $HOME/public_html/randomplots/saturation.pdf > out

# FDR calc
# comes from samocha set
bedtools intersect -a control.combine.vcf.gz -b essentials/gnomadbased-ccrs.bed.gz -wa -wb -sorted | awk '{print $NF}' > controlscores
bedtools intersect -a neurodev.combine.vcf.gz -b essentials/gnomadbased-ccrs.bed.gz -wa -wb -sorted | awk '{print $NF}' > neurodevscores
python fdr.py controlscores neurodevscores $HOME/public_html/randomplots/fdr.pdf

bedtools intersect -a <(zcat essentials/gnomadbased-ccrs.bed.gz | awk '$NF>=90') -b control.combine.vcf.gz -wa -wb -sorted | > benignbreaks.txt
vcfanno relative.conf control.combine.vcf.gz | bgzip -c > control.combine.ranges.vcf.gz
python relative.py control.combine.ranges.vcf.gz $HOME/public_html/randomplots/relative.pdf
