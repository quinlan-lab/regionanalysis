sed '1d' exacresiduals/results/2016_12_10/weightedresiduals.txt | cut -f -4,12 | sed 's/^/chr/g' > residuals.bed
liftOver residuals.bed $DATA/hg19ToDanRer10.over.chain.gz zebrafishregions.v10.bed unmapped #-positions argument prints coords in browser format
