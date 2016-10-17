# if you only want the top X number of genes from each list
numgenes=$1
# this created a gene rank for Wang et al based on an aggregate CRISPR score (CS) across all cell types (4)
python essentialgenes.py | sort -k2,2nr > intermeds/wanggenerank.txt
awk 'NR==FNR{a[$1]}; {for (i in a) if (i == $8) print}' <(head -$numgenes intermeds/wanggenerank.txt) FS="\"" <(zgrep -P 'protein_coding\tCDS|protein_coding\tstop_codon' $DATA/Homo_sapiens.GRCh37.75.gtf.gz) > genescreens/wanggenes.txt

# Shalem et al came pre-sorted by log2 mean depletion and only used one cell type KBM7 (log2=CS score in Wang et al)
awk 'NR==FNR{a[$1]}; {for (i in a) if (i == $8) print}' <(head -$numgenes ogfiles/shalemetal2014.txt) FS="\"" <(zgrep -P 'protein_coding\tCDS|protein_coding\tstop_codon' $DATA/Homo_sapiens.GRCh37.75.gtf.gz) > genescreens/shalemgenes.txt

# mouse KOs from Dickinson et al

tr -s "" "\n" < ogfiles/dickinsonetal2016.txt | grep -P "\tLethal" > genescreens/dickinson.txt

# intersection with residual percentiles

cut -f 1 genescreens/dickinson.txt | awk 'FNR==NR{a[$1];next}$4 in a{print $14}' - exacresiduals/results/weightedresiduals.txt > ryan/mousegenes

# pLI distro for top 0.1%

awk 'FNR==NR{a[$4];next}$2 in a{print $20}' <(awk '$14 > 99.9' exacresiduals/results/weightedresiduals.txt | sed '1d') $DATA/forweb_cleaned_exac_r03_march16_z_data_pLI.txt > ryan/pli

# to create wang.bed file for sgRNA analysis
tr -s "" "\n" < ogfiles/CSbysgRNA.txt | sed '1d' | cut -f -4 > genescreens/wang.bed 
