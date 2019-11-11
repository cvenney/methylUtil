#!/bin/bash

# get_TSS_and_promoters.sh

GFF="02_reference/genes.gff.gz"
GENOME="02_reference/genome.fasta"

if [[ ! -e 02_reference/genome.fasta.fai ]]
then
    samtools faidx 02_reference/genome.fasta
fi

cut -f1,2 02_reference/genome.fasta.fai > 02_reference/genome.genome

if [[ "${GFF##*.}" == "gz" ]]
then
	gunzip -c $GFF | awk '($1 !~ /#.*/ && $3 == "gene") {OFS="\t"} {if ($7 == "+") {print $1, $4 - 1, $4, $9, ".", $7} else if($7 == "-") {print $1, $5, $5 + 1, $9, ".", $7}}' > 05_bed_files_for_analysis/tss.bed
else
	cat $GFF | awk '($1 !~ /#.*/ && $3 == "gene") {OFS="\t"} {if ($7 == "+") {print $1, $4 - 1, $4, $9, ".", $7} else if($7 == "-") {print $1, $5, $5 + 1, $9, ".", $7}}' > 05_bed_files_for_analysis/tss.bed
fi

bedtools slop -i 05_bed_files_for_analysis/tss.bed -g $GENOME -l 1000 -r 200 -s > 05_bed_files_for_analysis/promotors.bed

