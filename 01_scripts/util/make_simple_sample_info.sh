#!/bin/bash
# script to generate a crude sample info file from the data

if [ $# -ne 1 ]
then
    echo "Usage: $0 <sample_info_prefix>"
    exit
fi

file=$1

if [ -e $file ]
then
    echo "File already exists!"
    exit
fi

echo "sample file" > $file


for i in $(ls 03_raw_bedGraphs/*.bedGraph.gz)
do
    sampleid=$(basename $i | perl -pe 's/\..*//')
    echo "${sampleid} 04_filtered_bedGraphs/$(basename ${i})" >> $file
done
