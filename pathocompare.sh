bedtools intersect -a <(sed '1d' exacresiduals/results/2016_11_17/weightedresiduals.txt) -b patho.vcf | cut -f 12 > ryan/patho
bedtools intersect -a <(sed '1d' exacresiduals/results/2016_11_17/weightedresiduals.txt) -b benign.vcf | cut -f 12 > ryan/benign
paste ryan/benign ryan/patho | python ryan/hist.py -o rh.png -t "Clinvar Variant Comparison"
