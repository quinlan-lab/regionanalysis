set -x
DATA=/scratch/ucgd/lustre/u1021864/serial
REGIONS=$1
HOTSPOT=$2
FOLDER=$3
OUT=$FOLDER/$FOLDER'_'$(echo $HOTSPOT | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")'.overlap.txt'

$HOME/software/poverlap/poverlap.py poverlap --a $REGIONS --b $HOTSPOT -g $DATA/human.hg19.genome --include $DATA/Homo_sapiens.GRCh37.75.gtf.gz -N 100 > results/$OUT
python plot.py results/$OUT $FOLDER
