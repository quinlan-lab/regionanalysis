DATA=/scratch/ucgd/lustre/u1021864/serial
HOME=/uufs/chpc.utah.edu/common/home/u1021864
#
python proplot.py -b exacresiduals/gnomad10x.5-ccrs.bed.gz -p $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170802.vcf.gz -c essentials/gnomad.vcf.gz -g exacresiduals/flatexome.bed.gz --region 20:62037542-62103993 --pfam pfam/pfam.exonic.bed.gz -f KCNQ2
bedtools intersect -a <(zcat exacresiduals/gnomad10x.5-ccrs.bed.gz | awk '$NF >= 90') -b <(awk '$NF <= 0.2' pli.bed) | bedtools intersect -a stdin -b <(zcat essentials/mpc.regions.clean.sorted.bed.gz | awk '$NF >= .6' ) | bedtools intersect -a stdin -b $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170802.vcf.gz > potentialb.bed
python proplot.py -b exacresiduals/gnomad10x.5-ccrs.bed.gz -p $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170802.vcf.gz -c essentials/gnomad.vcf.gz -g exacresiduals/flatexome.bed.gz --region 1:201328136-201346890 --pfam pfam/pfam.exonic.bed.gz -f TNNT2
