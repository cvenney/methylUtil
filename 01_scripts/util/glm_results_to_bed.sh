#!/bin/bash
# wald_results_to_bed.sh

if [ $# -ne 1 ]
then
    echo "Usage: $0 <glm_test_DMR_file>"
    exit
fi

file=$1
name=$(echo $file | perl -pe 's/(\.txt)(\.gz)//g')

if [ ${file##*.} == "gz" ]
then
    gunzip -c $file | 
    awk -v OFS="\t" '(NR != 1){
        print $1, $2, $3, "nCG=" $5 ";score=" $6, $6, "."
    }' | sort -k1,2 > ${name}.bed
else
    cat $file | 
    awk -v OFS="\t" '(NR != 1){
        print $1, $2, $3, "nCG=" $5 ";score=" $6, $6, "."
    }' | sort -k1,2 > ${name}.bed
fi

