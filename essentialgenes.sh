python essentialgenes.py > intermeds/wanggenerank.txt
awk 'NR==FNR{a[$1]}; {for (i in a) if (i == $8) print}' intermeds/wanggenerank.txt FS="\"" <(zcat $DATA/Homo_sapiens.GRCh37.75.sorted.gtf.gz | grep -e 'protein_coding\tCDS' -e 'protein_coding\tstop_codon') > wanggenes.txt
