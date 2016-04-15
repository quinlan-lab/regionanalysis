bedtools intersect -a - -b $DATA/vepcanonicalexons.gtf -wb | awk 'match($26, /\"(\w*)/, t) {if ($4 == t[1]) print}' | cut -f -10 | sort -k1,1 -k2,2n | uniq
