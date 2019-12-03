#!/bin/bash
# repeatmasker_to_bed.sh

if [ $# -ne 1 ]
then
    echo "Usage: $0 <repeatmasker_file.out>"
    exit
fi

outdir="05_bed_files"
file=$1
name=$(basename $(echo $file | perl -pe 's/\.out\.gz//'))

if [ ${file##*.} == "gz" ]
then
    gunzip -c $file | awk -v OFS="\t" '(NR > 3){
        print $5, $6, $7, $11 "=" $10
    }' > ${outdir}/${name}.bed
else
    cat $file | awk -v OFS="\t" '(NR > 3){
        print $5, $6, $7, $11 "=" $10
    }' > ${outdir}/${name}.bed
fi
