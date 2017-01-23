#get random genes for comparison (replace num with whatever number represents the top CCRs total; 837 for top 1%, 4476 for top 5%, 8399 for top 10%)
rm random.txt
file=$1
num=`wc -l $file | perl -pe 's/^\s*//g' | cut -f1 -d ' '`
max=100
for (( i=1; i <= $max; ++i ))
do
    cat genescreens/universe.tsv | ryan/reservoir_sampling $num > random_genes.txt
    echo "Wang et al	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 ryan/kbm7_score.txt | sort | uniq) | wc -l >> random.txt
    echo "ClinGen	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cat genescreens/clingen_level3_genes_2015_02_27.tsv | cut -f 1 | sort | uniq) | wc -l >> random.txt
    echo "Dickinson et al	randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 genescreens/dickinson.txt | sort | uniq) | wc -l >> random.txt
    #echo "Clinvar patho	randoms" >> random.txt
    #awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 genescreens/pathogenes.txt | sort | uniq) | wc -l >> random.txt
done
#get random averages
python randomparse.py random.txt > randomintersects.txt
