SITES=$1 REGIONS=$2 COMP=$3
HOME=/uufs/chpc.utah.edu/common/home/u1021864
DATA=/scratch/ucgd/lustre/u1021864/serial
set -e
check(){
        obs=$1
        exp=$2

if [[ `diff $obs $exp` != "" ]]; then
        echo "FAIL"
else
        echo "transcripts are correct"
fi
}
if [[ ! -f "$DATA/veptranscripts.txt" || ! -f "$DATA/vepcanonicalexons.gtf" ]]
then
    python vepcanon.py > $DATA/veptranscripts.txt
    cat <(zgrep '#!' $DATA/Homo_sapiens.GRCh37.75.gtf.gz) <(zgrep 'transcript_id' $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep -P 'protein_coding\tCDS|stop_codon\t') > foo.gtf
    awk -F "\"" 'NR==FNR {a[$1]} {for (i in a) if ($4 == i) print}' $DATA/veptranscripts.txt foo.gtf > $DATA/vepcanonicalexons.gtf
    check <(sort $DATA/veptranscripts.txt | uniq) <(grep -v ^# $DATA/vepcanonicalexons.gtf | cut -d "\"" -f 4 | sort | uniq)
fi
HOTSPOT=$HOME'/analysis/intermeds/'$(echo $COMP | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")'.gtf'
if [[ $SITES -eq 0 && ! -f "$HOTSPOT" ]]
then
    awk 'NR==FNR{a[$1]}; {for (i in a) if (i == $8) print}' $COMP FS="\"" $DATA/vepcanonicalexons.gtf > $HOTSPOT
fi
if [[ $SITES -eq 1 ]]
then
    HOTSPOT=$COMP
fi
FOLDER=$(echo $REGIONS | rev | cut -f 1 -d "/" | rev | perl -pe 's/(.*)gene|resid\..*/$1/g')
OUT=$FOLDER/$FOLDER'_'$(echo $HOTSPOT | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")'.overlap.txt'
if [[ ! -f "results/$OUT" ]]
then
    $HOME/software/poverlap/poverlap.py poverlap --a $REGIONS --b $HOTSPOT \
        -g $DATA/human.hg19.genome --include $DATA/vepcanonicalexons.gtf -N 1000 \
        > results/$OUT
else
    echo "results/$OUT already exists."
fi
