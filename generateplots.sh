set -x
DATA=/scratch/ucgd/lustre/u1021864/serial
REGIONS=regions/topresid.txt
# pass your file for montecarlo analysis, could be pathogenic variants, benign variants, whatever
HOTSPOT=$1
# first runs analysis on top residuals
FOLDER=top
function runplot {
  OUT=$FOLDER/$FOLDER'_'$(echo $HOTSPOT | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")'.overlap.txt'

  $HOME/software/poverlap/poverlap.py poverlap --a $REGIONS --b $HOTSPOT -g $DATA/human.hg19.genome --include $DATA/Homo_sapiens.GRCh37.75.gtf.gz -N 100 > results/$OUT
  python plot.py results/$OUT $FOLDER
}
# then runs analysis on expected regions, regions on the model line "middle residuals"
runplot
REGIONS=regions/midresid.txt
FOLDER=mid
runplot
