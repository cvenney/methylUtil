#!/bin/bash
# generate_genomic_context.sh

if [[ $# -ne 2 ]];then
	echo "Usage: $0 <\"WALD\" or \"GLM\"> <DML or DMR>"
	exit 1
fi

TEST=$1
INPUT=$2

if [[ $TEST != "WALD" ]] | [[ $TEST != "GLM" ]];then
    echo "ERROR: Test must be \"WALD\" or \"GLM\""
    exit 1
fi

if [[ $TEST == "WALD" ]];then
    01_scripts/util/wald_results_to_bed.sh $INPUT
fi

if [[ $TEST == "GLM" ]];then
    01_scripts/util/glm_results_to_bed.sh $INPUT
fi

BED="$(echo $INPUT | perl -pe 's/(\.txt)(\.gz)//g').bed"

GFF="02_reference/genes.gff.gz"

if [ ! -f $GFF ];then
    if [ ! -f ${GFF%.*} ];then
        echo "Error: GFF file not found!"
        exit
    else
    	GFF=${GFF%.*}
    fi
fi

if [ ! -e ${GFF%%.*}_with_utrs.gff ];then
	echo "Adding UTRs to GFF file..."
	python3 01_scripts/util/NCBI_add_utrs_to_gff.py $GFF > ${GFF%%.*}_with_utrs.gff
fi

echo "Intersecting DMRs with GFF..."
bedtools intersect -a $BED -b ${GFF%%.*}_with_utrs.gff -wao > ${BED%.*}_dmr_context.txt

echo "Printing GeneIDs..."
awk '($9 == "gene") {print $15}' ${BED%.*}_dmr_context.txt | \
	sed -e 's/ID=//' -e 's/;.*GeneID:/:/' -e 's/;.*//' |
	sort | uniq > ${BED%.*}_geneids.txt

echo "Summarizing results..."
01_scripts/util/get_genomic_context.R ${BED%.*}_dmr_context.txt
