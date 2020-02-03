#!/bin/bash
# split_by_chr.sh
# wrapped script to batch process bedGraphs in 02_data

if [ $# -ne 1 ]
then
    echo "Usage: $0 <chr_prefix (e.g. NC_036)>"
    exit
fi

prefix=$1

for file in $(ls 04_filtered_bedGraphs/*.bedGraph.gz)
do
    01_scripts/util/split_bedGraph_by_chr.py $file $prefix
done
