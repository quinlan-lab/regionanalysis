#mysql --user=genome --host=genome-mysql.soe.ucsc.edu -D hg19 -A -e  'SELECT chrom,chromStart,chromEnd,name from ucscGenePfam;' | perl -pe 's/^chr//' | bgzip -c > pfam.hg19.bed.gz
if [ ! -s pfamflat.bed ]; then
    python flattenpfams.py $DATA/pfam.bed > pfamflat.bed # flatten exomes across transcripts
fi
bedtools intersect -a ../exacresiduals/gnomad10x.5-ccrs.bed.gz -b pfamflat.bed -wao -sorted > ccr-pfam.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
cat ccr-pfam.bed | python ccrpfam.py
paste tmp/nodomccrs tmp/domccrs | python ../hist.py -o pfam_dist.pdf -t "gnomad based CCR v Pfam" -l "nodoms" "doms" # gnomad based CCR
awk '$14>90' ccr-pfam.bed | cut -f 18 | sort | uniq -c | sort -k1,1nr | head -100 | sed '1d' | sed 's/^\s*\S*//g' > top100doms #top 100 domains
bedtools intersect -b ../exacresiduals/gnomad10x.5-ccrs.bed.gz -a pfamflat.bed -u -sorted | cut -f 4 | sort | uniq > pfams.txt
#5117/5416 uniq pfams covered by gnomad-ccrs.bed.gz
bedtools intersect -a pfamflat.bed -b ../exacresiduals/gnomad10x.5-ccrs.bed.gz -wb -sorted | sort -k4,4 > pfam-ccr.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
python fameval.py -i pfam-ccr.bed > pfamshist.txt
python plotpfam.py -p top100doms -s pfamshist.txt -o $HOME/public_html/randomplots/pfam_hists.pdf
#python plotpfam.py -p pfams.txt -s pfamshist.txt -o pfam_hists.pdf
python gerppfam.py
