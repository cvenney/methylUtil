#!/bin/bash

# remove_CT_SNPs_from_bedGraphs.sh

INPUT="03_methylation_bedGraphs"
OUTPUT="04_filtered_bedGraphs"
NCPUS=4

# Negative intersect with bedtools
for file in $(ls "$INPUT"/*.bedGraph.gz | perl -pe 's/\.bedGraph\.gz//g')
do
    name=$(basename $file)
    
    echo "Filtering sample: $file"    
    
    gunzip -c "$INPUT"/"$name".bedGraph.gz | 
    bedtools intersect -a stdin -b 02_reference/CT_AG_snps.bed -wa -v | 
    gzip > "$OUTPUT"/"$name".bedGraph.gz

done
