#!/bin/bash
REGIONS=$1
for i in ogfiles/*; do
    if [[ $(awk 'END {print NF}' $i) -le 2 ]]
    then
        bash analyze.sh 0 $REGIONS $i
    else
        bash analyze.sh 1 $REGIONS $i
    fi
done 
