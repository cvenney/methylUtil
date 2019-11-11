#!/bin/bash

# remove_CT_SNPs_from_bedGraphs.sh

INPUT="03_methylation_bedGraphs"
OUTPUT="04_filtered_bedGraphs"
NCPUS=4

# Negative intersect with bedtools
ls "$INPUT"/*.bedGraph.gz | perl -pe 's/\.bedGraph\.gz//g' |
parallel --bar -j $NCPUS gunzip -c {}.bedGraph.gz | 
bedtools intersect -a stdin -b 02_reference/CT_AG_snps.bed -wa -v | 
gzip > "$OUTPUT"/{}.bedGraph.gz

