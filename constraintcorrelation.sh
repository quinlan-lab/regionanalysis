# intersect CCRs >=95 with pLI, missense Z and RVIS, match on gene
bedtools intersect -a $HOME/software/pathoscore/score-sets/GRCh37/mis_Z/missensez.bed.gz -b <(zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5 >= 95') -wao -sorted > miszccr.txt
bedtools intersect -a $HOME/software/pathoscore/score-sets/GRCh37/RVIS/rvis.bed.gz -b <(zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5 >= 95') -wao -sorted > rvisccr.txt
bedtools intersect -a $HOME/software/pathoscore/score-sets/GRCh37/pLI/pLI.bed.gz -b <(zcat $HOME/public_html/files/ccrs.v2.20180420.bed12.bed.gz | awk '$5 >= 95') -wao -sorted > pliccr.txt
# create correlation calc and scatter plots
python conscorr.py miszccr.txt rvisccr.txt pliccr.txt
