top=$1
date=2016_12_10
#remove previous files 
rm intersects.txt
#create top 1% gene list, 5%, 10% etc.
cut -f 4 <(awk -v t=$top 'BEGIN {top=100-t} {if ($12>top) print}' exacresiduals/results/$date/weightedresiduals.txt) | sed '1d' | sort | uniq | awk '{split($1,a,","); for (i in a) print a[i];}' > /tmp/resids
wc -l /tmp/resids
#compare top gene list to an essential genescreen file of some kind
echo "ClinGen	total" >> intersects.txt
wc -l genescreens/clingen_level3_genes_2015_02_27.tsv >> intersects.txt
echo "ClinGen	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' /tmp/resids <(cat genescreens/clingen_level3_genes_2015_02_27.tsv | cut -f 1 | sort | uniq) | wc -l >> intersects.txt
echo "Dickinson et al	total" >> intersects.txt
wc -l genescreens/dickinson.txt >> intersects.txt
echo "Dickinson et al	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' /tmp/resids <(cat genescreens/dickinson.txt | cut -f 1 | sort | uniq) | wc -l >> intersects.txt
echo "Wang et al	total" >> intersects.txt
cut -f 1 ryan/kbm7_score.txt | sort | uniq | wc -l >> intersects.txt
echo "Wang et al	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' /tmp/resids <(cut -f 1 ryan/kbm7_score.txt | sort | uniq) | wc -l >> intersects.txt
#getting clinvar pathogenic intersections...use pathogenes.py on the clinvar pathogenic vcf
#python pathogenes.py patho.vcf
#echo "Clinvar pathogenic	total" >> intersects.txt
#wc -l genescreens/pathogenes.txt >> intersects.txt
#echo "Clinvar pathogenic	intersections" >> intersects.txt
#awk 'FNR==NR{a[$1]; next}$1 in a{print}' resids <(cut -f 1 genescreens/pathogenes.txt | sort | uniq) | wc -l >> intersects.txt
#python genebar.py genefractions
