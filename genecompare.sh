top=$1
#remove previous files 
rm intersects.txt; rm random.txt
#create top 1% gene list, 5%, 10% etc.
cut -f 4 <(awk -v t=$top 'BEGIN {top=100-t} {if ($12>top) print}' exacresiduals/results/2016_11_17/nosizefilter/weightedresiduals.txt) | sed '1d' | sort | uniq | awk '{split($1,a,","); for (i in a) print a[i];}' > foo
#compare top gene list to an essential genescreen file of some kind
echo "pLI>.9	total" >> intersects.txt
awk '$20>.9' $DATA/*pLI* | cut -f 2 | sed '1d' | sort | uniq | wc -l >> intersects.txt
echo "pLI>.9	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' foo <(awk '$20>.9' $DATA/*pLI* | cut -f 2 | sed '1d' | sort | uniq) | wc -l >> intersects.txt
echo "Clingen	total" >> intersects.txt
wc -l genescreens/clingen_level3_genes_2015_02_27.tsv >> intersects.txt
echo "Clingen	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' foo <(cat genescreens/clingen_level3_genes_2015_02_27.tsv | cut -f 1 | sort | uniq) | wc -l >> intersects.txt
echo "Dickinson	total" >> intersects.txt
wc -l genescreens/dickinson.txt >> intersects.txt
echo "Dickinson	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' foo <(cat genescreens/dickinson.txt | cut -f 1 | sort | uniq) | wc -l >> intersects.txt
echo "Wang	total" >> intersects.txt
cut -f 1 ryan/kbm7_score.txt | sort | uniq | wc -l >> intersects.txt
echo "Wang	intersections" >> intersects.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' foo <(cut -f 1 ryan/kbm7_score.txt | sort | uniq) | wc -l >> intersects.txt
#getting clinvar pathogenic intersections...use pathogenes.py on the clinvar pathogenic vcf
#python pathogenes.py patho.vcf
#echo "Clinvar pathogenic	total" >> intersects.txt
#wc -l genescreens/pathogenes.txt >> intersects.txt
#echo "Clinvar pathogenic	intersections" >> intersects.txt
#awk 'FNR==NR{a[$1]; next}$1 in a{print}' foo <(cut -f 1 genescreens/pathogenes.txt | sort | uniq) | wc -l >> intersects.txt
#get random genes for comparison (replace 837 with whatever number represents the top CCRs total; 837 for top 1%, 4476 for top 5%, 8399 for top 10%)
num=`wc -l foo | cut -f1 -d' '`
max=100
for (( i=1; i <= $max; ++i ))
do
    cat genescreens/universe.tsv | ryan/reservoir_sampling $num > random_genes.txt
    echo "pLI>.9	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(awk '$20>.9' $DATA/*pLI* | cut -f 2 | sed '1d' | sort | uniq) | wc -l >> random.txt
    echo "Dickinson	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 genescreens/dickinson.txt | sort | uniq) | wc -l >> random.txt
    echo "Wang	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 ryan/kbm7_score.txt | sort | uniq) | wc -l >> random.txt
    echo "Clingen	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cat genescreens/clingen_level3_genes_2015_02_27.tsv | cut -f 1 | sort | uniq) | wc -l >> random.txt
    #echo "Clinvar patho	randoms" >> random.txt
    #awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 genescreens/pathogenes.txt | sort | uniq) | wc -l >> random.txt
done
#get random averages
python randomparse.py random.txt > randomintersects.txt
#python genebar.py genefractions
