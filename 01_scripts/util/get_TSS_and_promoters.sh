#!/bin/bash

# get_TSS_and_promoters.sh

GFF="02_reference/genes.gff.gz"
GENOME="02_reference/genome.genome"

if [[ ! -e 02_reference/genome.fasta.fai ]]
then
    samtools faidx 02_reference/genome.fasta
fi

cut -f1,2 02_reference/genome.fasta.fai > 02_reference/genome.genome

if [[ "${GFF##*.}" == "gz" ]]
then
	gunzip -c $GFF | 
	awk -v OFS="\t" '($1 !~ /"#.*"/ && $3 == "gene") {
	    if ($7 == "+") {
	        print $1, $4 - 1, $4, $9, ".", $7
	    } else if($7 == "-") {
	        print $1, $5, $5 + 1, $9, ".", $7}}' | 
	perl -pe 's/ID=//' | 
	perl -pe 's/;Dbxref=GeneID//' | 
	perl -pe 's/;.*\b//' > 05_bed_files/tss.bed
else
	cat $GFF | 
	awk -v OFS="\t" '($1 !~ /"#.*"/ && $3 == "gene") {
	    if ($7 == "+") {
	        print $1, $4 - 1, $4, $9, ".", $7
	    } else if($7 == "-") {
	        print $1, $5, $5 + 1, $9, ".", $7}}' | 
	perl -pe 's/ID=//' | 
	perl -pe 's/;Dbxref=GeneID//' | 
	perl -pe 's/;.*\b//' > 05_bed_files/tss.bed
fi



bedtools slop -i 05_bed_files/tss.bed -g $GENOME -l 1000 -r 200 -s > 05_bed_files/promotors.bed

