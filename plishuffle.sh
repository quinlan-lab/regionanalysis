rm random.txt
num=`wc -l foo | cut -f1 -d' '`
max=100
for (( i=1; i <= $max; ++i ))
do
    cat genescreens/universe.tsv | ryan/reservoir_sampling $num > random_genes.txt
    echo "pLI > .9 randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(awk '$20>.9' $DATA/*pLI* | cut -f 2 | sed '1d' | sort | uniq) | wc -l >> random.txt
    echo "Dickinson randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 genescreens/dickinson.txt | sort | uniq) | wc -l >> random.txt
    echo "Wang randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 ryan/kbm7_score.txt | sort | uniq) | wc -l >> random.txt
    echo "Clingen randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cat genescreens/clingen_level3_genes_2015_02_27.tsv | cut -f 1 | sort | uniq) | wc -l >> random.txt
    echo "Clinvar patho randoms" >> random.txt
    awk 'FNR==NR{a[$1]; next}$1 in a{print}' random_genes.txt <(cut -f 1 genescreens/pathogenes.txt | sort | uniq) | wc -l >> random.txt
done
