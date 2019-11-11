#!/bin/bash

# find_CpG_islands.sh

cpgplot -sequence 02_reference/genome.fasta \
    -window 100 \
    -minlen 200 \
    -minoe 0.6 \
    -minpc 50. \
    -outfile 05_bed_files_for_analysis/genome.cpgplot \
    -noplot \
    -outfeat 05_bed_files_for_analysis/cpg_islands.gff
    
if [[ ! -e 02_reference/genome.fasta.fai ]]
then
    samtools faidx 02_reference/genome.fasta
fi

cut -f1,2 02_reference/genome.fasta.fai > 02_reference/genome.genome

bedtools merge -i 05_bed_files_for_analysis/cpg_islands.gff > 05_bed_files_for_analysis/cpg_islands.bed
rm 05_bed_files_for_analysis/cpg_islands.gff

bedtools slop -i 05_bed_files_for_analysis/cpg_islands.bed -g 02_reference/genome.genome -b 2000 |
bedtools merge -i stdin |
bedtools subtract -a stdin -b 05_bed_files_for_analysis/cpg_islands.bed > 05_bed_files_for_analysis/cpg_shores.bed

