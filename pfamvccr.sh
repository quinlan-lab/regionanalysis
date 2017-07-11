if [ ! -s pfamflat.bed ]; then
    python flattenpfams.py $DATA/pfam.bed > pfamflat.bed # flatten exomes across transcripts
fi
bedtools intersect -a gnomad-ccrs.bed.gz -b pfamflat.bed -wao -sorted > foo # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
cat foo | python ccrpfam.py
paste tmp/nodomccrs tmp/domccrs | python hist.py -o pfam_dist.pdf -t "gnomad based CCR v Pfam" -l "nodoms" "doms" # gnomad based CCR
