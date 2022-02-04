#!/bin/bash

# remove_CT_SNPs_from_bedGraphs.sh

INPUT="/project/lbernatchez/01_projects/Projects/Clare_Venney/01_bwa-meth_pipeline_brook_charr_epigenetics/07_methyl_dackel"
OUTPUT="04_filtered_bedGraphs"
NCPUS=4

# Negative intersect with bedtools
for file in $(ls "$INPUT"/*.bedGraph.gz | perl -pe 's/\.bedGraph\.gz//g')
do
    name=$(basename $file)
    
    echo "Filtering sample: $file"    
    
    gunzip -c "$INPUT"/"$name".bedGraph.gz | 
    bedtools intersect -a stdin -b 02_reference/HI.5144.008.IDT_i7_86---IDT_i5_86.Sfon_WGS.dedup.no_overlap_freebayes_CT_AG_snps.bed -wa -v | 
    gzip > "$OUTPUT"/"$name".bedGraph.gz

done
