#!/bin/bash

# remove_CT_SNPs_from_bedGraphs.sh

INPUT="/project/lbernatchez/01_projects/Projects/Clare_Venney/01_bwa-meth_pipeline_brook_charr_epigenetics/07_methyl_dackel"
OUTPUT="04_filtered_bedGraphs"
NCPUS=4

# Negative intersect with bedtools
for file in $(ls "$INPUT"/*.methylKit.gz | perl -pe 's/\.methylKit\.gz//g')
do
    name=$(basename $file)
    
    echo "Filtering sample: $file"    
    
    gunzip -c "$INPUT"/"$name".methylKit.gz | 
    awk 'BEGIN{OFS="\t"}(NR!=1){print $2, $3 -1, $3, $1, $4, $5, $6, $7}' |
    bedtools intersect -a stdin -b 02_reference/HI.5144.008.IDT_i7_86---IDT_i5_86.Sfon_WGS.dedup.no_overlap_freebayes_CT_AG_snps.bed -wa -v |
    awk 'BEGIN{OFS="\t"; print "chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"} {print $4, $1, $3, $5, $6, $7, $8}' |
    gzip > "$OUTPUT"/"$name".methylKit.gz

done
