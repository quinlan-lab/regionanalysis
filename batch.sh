#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o %j-%N.out
#SBATCH -e %j-%N.err
#SBATCH --time=4:00:00
FILE=$1;FOLDER=$2
python cpgexoncalc.py
#bash loop.sh $FILE
#bash plotloop.sh $FOLDER
#bash analyze.sh 0 regions/bottomgerpresid.txt ogfiles/all_ar.tsv
#python plot.py results/all_ar.overlap.txt
#python regression.py
#python plotdistro.py
