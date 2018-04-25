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
HOME=/uufs/chpc.utah.edu/common/home/u1021864
DATA=/scratch/ucgd/lustre/u1021864/serial
if [ ! -s pfam.exonic.bed ]; then
    sed 's/^chr//g' pfam.browser.bed | awk '{split($4,a,"_exon"); print $1, $2, $3, a[1]}' OFS='\t' | grep -P "^1|^2|^3|^4|^5|^6|^7|^8|^9|^10|^11|^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22"| bedtools intersect -a stdin -b $HOME/analysis/exacresiduals/flatexome.bed -wb | awk '{print $1,$2,$3,$4,$8}' OFS='\t' | sort -k4,4 -k1,1 -k2,2n > pfam.exonic.bed
    sort -k1,1 -k2,2n pfam.exonic.bed | bgzip -c > pfam.exonic.bed.gz; tabix pfam.exonic.bed.gz
    #python flattenpfams.py $DATA/pfam.bed > pfam.exonic.bed # flatten exomes across transcripts
fi
if [ ! -s pfam.genome.gene.bed ]; then
    sed 's/^chr//g' pfam.genome.bed | grep -P "^1|^2|^3|^4|^5|^6|^7|^8|^9|^10|^11|^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22"| bedtools intersect -a stdin -b $HOME/analysis/exacresiduals/flatexome.bed -wa -wb | awk '{print $1,$2,$3,$4,$8}' OFS='\t' | sort -k4,4 -k1,1 -k2,2n | uniq > pfam.genome.gene.bed
    sort -k1,1 -k2,2n pfam.genome.gene.bed | bgzip -c > pfam.genome.gene.bed.gz; tabix pfam.genome.gene.bed.gz
    #python flattenpfams.py $DATA/pfam.bed > pfam.exonic.bed # flatten exomes across transcripts
fi
if [ ! -s Pfam-A.clans.tsv ]; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
    gunzip Pfam-A.clans.tsv.gz
fi
bedtools intersect -a ../exacresiduals/gnomad10x.5syn-ccrs.bed.gz -b pfam.genome.gene.bed -wao > ccr-pfam.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
cat ccr-pfam.bed | python ccrpfam.py
paste tmp/nodomccrs tmp/domccrs | python ../hist.py -o pfam_dist.pdf -t "gnomad based CCR v Pfam" -l "nodoms" "doms" # gnomad based CCR
awk '$14>90' ccr-pfam.bed | cut -f 18 | sort | uniq -c | sort -k1,1nr | head -200 | sed '1d' | sed 's/^\s*\S*//g' > top100doms #top 100 domains
bedtools intersect -b ../exacresiduals/gnomad10x.5syn-ccrs.bed.gz -a pfam.genome.gene.bed -u | cut -f 4 | sort | uniq > pfams.txt # all ccr intersecting pfams
cut -f 4 pfam.genome.gene.bed | sort | uniq -c | sed 's/^\s*//g' > pfamcounts.txt # counts of pfams
#
# histograms of ccr dists across pfams #
#
bedtools intersect -a pfam.genome.gene.bed -b ../exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb | sort -k4,4 > pfam-ccr.bed # pfam.bed is GRCh37.gtf with intron-containing pfam file across chrom 1-22, X, Y;
python fameval.py -i pfam-ccr.bed > pfamshist.txt
python plotpfam.py -p pfams.txt -s pfamshist.txt -c Pfam-A.clans.tsv -q pfamcounts.txt -o $HOME/public_html/randomplots/pfam_hists\(supp_doc_1\).pdf
#
# generate gerp v ccr plots
# 
sort -k4,4 pfam.genome.bed > pfamsorted.bed
python gerppfam.py pfamsorted.bed $DATA/hg19.gerp.bw $HOME/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz
python plotgerp.py ccrgerppfam.pkl $HOME/public_html/randomplots/gerpvccr_pfam.pdf
python ccrvgerp.py $DATA/hg19.gerp.bw $HOME/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz
python plotgerp.py ccrgerp.pkl $HOME/public_html/randomplots/gerpvccr.pdf purifyingselectionregions\(supp_table_2\).tsv
#
# nodom analysis
#
zcat $HOME/analysis/essentials/gnomadbased-ccrs.bed.gz | awk '{n=split($7,a,/,/); split(a[1],s,/-/); m=split(a[n],e,/-/); {printf "%s\t%s\t%s", $1, s[1], e[m]} {for (x=4; x<=NF; x++) printf "\t%s", $x} {printf "\n"}}' OFS="\t" | uniq > ccrsingenomespace.bed
bedtools intersect -a <(awk '$NF >= 95' ccrsingenomespace.bed) -b pfam.genome.gene.bed -v > topnonpfamccrs.bed
bedtools intersect -a <(awk '$NF <= 20 && $NF != 0' ccrsingenomespace.bed) -b pfam.genome.gene.bed -v > bottomnonpfamccrs.bed
bedtools intersect -a <(awk '$NF >= 95' ccrsingenomespace.bed) -b pfam.genome.gene.bed -u > toppfamccrs.bed
python nodomplot.py topnonpfamccrs.bed bottomnonpfamccrs.bed pfam.genome.gene.bed $HOME/analysis/exacresiduals/flatexome.bed $HOME/public_html/randomplots/nodom.pdf > output

# generate table of:
TOT=$(wc -l toppfamccrs.bed)
# 1. CCRs that enclose a Pfam domain
CP=$(bedtools intersect -a toppfamccrs.bed -b pfam.genome.bed -wa -wb | awk '{if ($3 >= $17 && $2 <= $16) print}' OFS="\t" | cut -f 1,2,3 | uniq | wc -l)
# 2. Pfam domains that enclose a CCR
PC=$(bedtools intersect -a toppfamccrs.bed -b pfam.genome.bed -wa -wb | awk '{if ($3 <= $17 && $2 >= $16) print}' OFS="\t" | cut -f 1,2,3 | uniq | wc -l)
# 3. CCRs that partially overlap a Pfam domain
PO=$(bedtools intersect -a toppfamccrs.bed -b pfam.genome.bed -wa -wb | awk '{if (($3 <= $17 && $2 <= $16) || ($3 >= $17 && $2 >= $16)) print}' OFS="\t" | cut -f 1,2,3 | uniq | wc -l)
bc <<< "$CP+$PC+$PO"
echo $TOT
# 4. CCRs that are near a Pfam domain  (need to define breakpoint)
# 5. CCRs that are far away from a Pfam domain (need to define breakpoint)
# from nodomplot.py:

# high ccrs not intersecting with PFam domains
bedtools intersect -a <(zcat $HOME/public_html/files/ccrs.v1.20171112.bed12.bed.gz | awk '$5>=90') -b pfam.genome.gene.bed.gz -sorted -v | wc -l
bedtools intersect -a <(zcat $HOME/public_html/files/ccrs.v1.20171112.bed12.bed.gz | awk '$5>=99') -b pfam.genome.gene.bed.gz -sorted -v | wc -l
