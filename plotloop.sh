#!/bin/bash
FOLDER=$1
for i in results/$FOLDER/*; do
    python plot.py $i $FOLDER
done
