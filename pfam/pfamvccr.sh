#mysql --user=genome --host=genome-mysql.soe.ucsc.edu -D hg19 -A -e  'SELECT chrom,chromStart,chromEnd,name from ucscGenePfam;' | perl -pe 's/^chr//' | bgzip -c > pfam.hg19.bed.gz
##############################################################################################
# make input: "Pfam A domains with chromosomal co-ordinates from UCSC table browser" #pfam.browser.bed
# Dated: 17April 2013
# -> Go to UCSC table browser:
# -> choose; assembly: hg19, genome: Human, region: genome
# -> Select 'Genes and Gene Prediction Tracks' under group and 'Pfam in UCSC Gene' under track
# -> select 'BED- browser extensible data' under output format and click 'get output'
# -> In the resulting page, select 'Exons plus' under 'Create one BED record per:' and 
# -> click get BED
##############################################################################################
if [ ! -s pfam.exonic.bed ]; then
    sed 's/^chr//g' pfam.browser.bed | awk '{split($4,a,"_exon"); print $1, $2, $3, a[1]}' OFS='\t' | grep -P "^1|^2|^3|^4|^5|^6|^7|^8|^9|^10|^11|^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22"| bedtools intersect -a stdin -b $HOME/analysis/exacresiduals/flatexome.bed -wb | awk '{print $1,$2,$3,$4,$8}' OFS='\t' | sort -k4,4 -k1,1 -k2,2n > pfam.exonic.bed
    sort -k1,1 -k2,2n pfam.exonic.bed | bgzip -c > pfam.exonic.bed.gz; tabix pfam.exonic.bed.gz
    #python flattenpfams.py $DATA/pfam.bed > pfam.exonic.bed # flatten exomes across transcripts
fi
if [ ! -s Pfam-A.clans.tsv ]; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
    gunzip Pfam-A.clans.tsv.gz
fi
bedtools intersect -a ../exacresiduals/gnomad10x.5-ccrs.bed.gz -b pfam.exonic.bed -wao > ccr-pfam.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
cat ccr-pfam.bed | python ccrpfam.py
paste tmp/nodomccrs tmp/domccrs | python ../hist.py -o pfam_dist.pdf -t "gnomad based CCR v Pfam" -l "nodoms" "doms" # gnomad based CCR
awk '$14>90' ccr-pfam.bed | cut -f 18 | sort | uniq -c | sort -k1,1nr | head -100 | sed '1d' | sed 's/^\s*\S*//g' > top100doms #top 100 domains
bedtools intersect -b ../exacresiduals/gnomad10x.5-ccrs.bed.gz -a pfam.exonic.bed -u | cut -f 4 | sort | uniq > pfams.txt
#5117/5416 uniq pfams covered by gnomad-ccrs.bed.gz
bedtools intersect -a pfam.exonic.bed -b ../exacresiduals/gnomad10x.5-ccrs.bed.gz -wb | sort -k4,4 > pfam-ccr.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
python fameval.py -i pfam-ccr.bed > pfamshist.txt
python plotpfam.py -p top100doms -s pfamshist.txt -c Pfam-A.clans.tsv -o $HOME/public_html/randomplots/pfam_hists.pdf
#python plotpfam.py -p pfams.txt -s pfamshist.txt -o pfam_hists.pdf
python gerppfam.py
python plotgerp.py ccrgerppfam.pkl
python ccrvgerp.py
python plotgerp.py ccrgerp.pkl
bedtools intersect -a $HOME/analysis/essentials/gnomadbased-ccrs.bed.gz -b pfam.exonic.bed.gz -v | awk '$14 >= 95' > topnonpfamccrs.bed
python nodomplot.py topnonpfamccrs.bed pfam.exonic.bed.gz $HOME/analysis/exacresiduals/flatexome.bed $HOME/public_html/randomplots/nodom.pdf > output
