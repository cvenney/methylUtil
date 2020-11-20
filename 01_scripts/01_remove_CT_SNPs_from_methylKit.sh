#!/bin/bash

# remove_CT_SNPs_from_bedGraphs.sh

INPUT="03_raw_bedGraphs"
OUTPUT="04_filtered_bedGraphs"
NCPUS=4

# Negative intersect with bedtools
for file in $(ls "$INPUT"/*.methylKit.gz | perl -pe 's/\.methylKit\.gz//g')
do
    name=$(basename $file)
    
    echo "Filtering sample: $file"    
    
    gunzip -c "$INPUT"/"$name".methylKit.gz | 
    awk 'BEGIN{OFS="\t"}(NR!=1){print $2, $3 -1, $3, $1, $4, $5, $6, $7}' |
    bedtools intersect -a stdin -b 02_reference/CT_AG_snps.bed -wa -v |
    awk 'BEGIN{OFS="\t"; print "chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"} {print $4, $1, $3, $5, $6, $7, $8}' |
    gzip > "$OUTPUT"/"$name".methylKit.gz

done
