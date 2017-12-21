#wget https://raw.githubusercontent.com/recluze/deepseq/master/predictions-unknowns.txt
#wget http://fastademo.bioch.virginia.edu/pfam_dna/*
#cp ~/analysis/essentials/pfam.bed .
awk 'FNR==NR{a[$1]=$2; next} {split($19,b,"\""); for (i in a) if (i==b[2]) print $0, a[i]}' FS=" " predictions-unknowns.txt pfam.bed > unknowns.bed
bedtools intersect -a ../../exacresiduals/gnomad10x.5-ccrs.bed.gz -b unknowns.bed -sorted -wb | sort -k14,14nr > ccrs-unknowns.bed
python pfamhist.py ccrs-unknowns.bed
