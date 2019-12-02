#!/bin/bash
# script to generate a crude sample info file from the data

if [ $# -ne 1 ]
then
    echo "Usage: $0 <sample_info_prefix>"
fi

file=$1

if [ -e $file ]
then
    echo "File already exists!"
    exit
fi

echo "sample\tfile" > $file


for i in $(ls 03_raw_bedGraphs/*.bedGraph.gz)
do
    sampleid=$(echo $file | perl -pe 's/*Index_[0-9]+\.//' | perl -pe 's/\.*//')
    echo "${sampleid}\t${i}" >> $file
done
