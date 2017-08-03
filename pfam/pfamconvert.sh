sed 's/^chr//g' pfam.browser.bed | awk '{split($4,a,"_exon"); print $1, $2, $3, a[1]}' OFS='\t' | sort -k4.4 -k1,1 -k2,2n > pfam.exonic.bed
