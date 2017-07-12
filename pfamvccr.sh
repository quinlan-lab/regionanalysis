if [ ! -s pfamflat.bed ]; then
    python flattenpfams.py $DATA/pfam.bed > pfamflat.bed # flatten exomes across transcripts
fi
bedtools intersect -a gnomad-ccrs.bed.gz -b pfamflat.bed -wao -sorted > ccr-pfam.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
cat ccr-pfam.bed | python ccrpfam.py
paste tmp/nodomccrs tmp/domccrs | python hist.py -o pfam_dist.pdf -t "gnomad based CCR v Pfam" -l "nodoms" "doms" # gnomad based CCR
awk '$14>90' ccr-pfam.bed | cut -f 18 | sort | uniq -c | sort -k1,1nr | head -11 | sed '1d' | sed 's/^\s*\S*//g' > top10doms #top 10 domains
bedtools intersect -a pfamflat.bed -b gnomad-ccrs.bed.gz -wb -sorted | sort -k4,4 > pfam-ccr.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
python fameval.py -i pfam-ccr.bed > pfamshist.txt
python plotpfam.py -p top10doms -b pfamshist.txt -o pfam_hists.pdf

