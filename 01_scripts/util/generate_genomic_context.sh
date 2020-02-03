#!/bin/bash
# generate_genomic_context.sh

if [ $# -ne 1 ];then
	echo "Usage: $0 <DML or DMR BED file>"
fi

BED=$1

GFF="02_reference/genes.gff.gz"

if [ ! -e ${GFF%%.*}_with_utrs.gff ];then
	echo "Adding UTRs to GFF file..."
	python3 01_scripts/util/NCBI_add_utrs_to_gff.py $GFF > ${GFF%%.*}_with_utrs.gff
fi

echo "Intersecting DMRs with GFF..."
bedtools intersect -a $BED -b ${GFF%%.*}_with_utrs.gff -wao > ${BED%.*}_dmr_context.txt

echo "Summarizing results..."
01_scripts/util/get_genomic_context.R ${BED%.*}_dmr_context.txt
