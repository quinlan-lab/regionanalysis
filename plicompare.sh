# for pLI comparison to essential gene lists
rm pli.txt
awk '$20>.9' $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt | cut -f 2 | sed '1d' | sort | uniq > /tmp/pli
echo "ClinGen	pLI" >> pli.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' /tmp/pli <(cat genescreens/clingen_level3_genes_2015_02_27.tsv | cut -f 1 | sort | uniq) | wc -l >> pli.txt
echo "Dickinson et al	pLI" >> pli.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' /tmp/pli <(cat genescreens/dickinson.txt | cut -f 1 | sort | uniq) | wc -l >> pli.txt
echo "Wang et al	pLI" >> pli.txt
awk 'FNR==NR{a[$1]; next}$1 in a{print}' /tmp/pli <(cut -f 1 ryan/kbm7_score.txt | sort | uniq) | wc -l >> pli.txt
