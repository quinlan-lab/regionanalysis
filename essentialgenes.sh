numgenes=$1
#python essentialgenes.py | sort -k2,2nr > intermeds/wanggenerank.txt
awk 'NR==FNR{a[$1]}; {for (i in a) if (i == $8) print}' <(head -$numgenes intermeds/wanggenerank.txt) FS="\"" <(zgrep -P 'protein_coding\tCDS|protein_coding\tstop_codon' $DATA/Homo_sapiens.GRCh37.75.gtf.gz) > genescreens/wanggenes.txt
