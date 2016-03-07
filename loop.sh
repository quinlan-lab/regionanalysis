#!/bin/bash
cd regions
FILES=$(ls)
for f in $FILES
do
    for g in $FILES
    do
        python ~/analysis/maketable.py $f $g 1 
    done
done
