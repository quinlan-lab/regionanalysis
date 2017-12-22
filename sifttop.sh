# filter out simple repeats and segmental duplications
bedtools intersect -v -a <(sed '1d' exacresiduals/fullresiduals.txt | bash dups.sh | awk '$12 > 99') -b $DATA/hgsegmental.bed | bedtools intersect -a stdin -b $DATA/hgsimple.bed -sorted -wao | awk '{t=$3-$2+1} {if ($NF/t < .2) print}' | cut -f -13 | sort -k1,1 -k2,2n | uniq > toprfiltered.txt
#look to see if there are non-PASS variants
cat <(zgrep '^#' $DATA/ExAC.r1.sites.vep.tidy.vcf.gz) <(zgrep -v 'PASS' $DATA/ExAC.r0.3.sites.vep.tidy.vcf.gz | grep -v '^#' | sort -k1,1 -k2,2n) > /tmp/f1l3.vcf
python nonpass.py /tmp/f1l3.vcf > nonpass.vcf
bedtools intersect -a toprfiltered.txt -b nonpass.vcf -sorted -wao | awk '{t=$3-$2+1} {for(i=1;i<=13;i++){printf "%s\t", $i}} {print $NF/t}' OFS='\t' | bedtools intersect -a stdin -b patho.vcf -v -sorted | uniq > nonpasstop.txt
#to check any protein use REST API: e.g., http://www.rcsb.org/pdb/protein/LMNA?chromosome=1&position=156105817-156106009
