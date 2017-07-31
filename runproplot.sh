DATA=/scratch/ucgd/lustre/u1021864/serial
HOME=/uufs/chpc.utah.edu/common/home/u1021864
python proplot.py -t ENST00000359125 -b exacresiduals/gnomad10x.5-ccrs.bed.gz -p $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -c essentials/gnomad.vcf.gz -g $DATA/Homo_sapiens.GRCh37.75.gtf.gz --region 20:62037542-62103993 --pfam $DATA/pfam.bed.gz -f KCNQ2
python proplot.py -t ENST00000245503 -b exacresiduals/gnomad10x.5-ccrs.bed.gz -p $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz -c essentials/gnomad.vcf.gz -g $DATA/Homo_sapiens.GRCh37.75.gtf.gz --region 17:10424600-10451237 --pfam $DATA/pfam.bed.gz -f MYH2
