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

GFF="02_reference/genes.gff"

if [ ! -f $GFF ];then
    echo "Error: GFF file not found!"
    exit
fi


echo "Intersecting DMRs with GFF..."
bedtools window -a ${GFF} -b $BED -l 5000 -r 5000 -sw > ${BED%.*}_dmr_context.txt

echo "Printing GeneIDs..."
awk '($3 == "gene") {print $9}' ${BED%.*}_dmr_context.txt | \
	sed -e 's/ID=//' -e 's/;.*GeneID:/:/' -e 's/;.*//' |
	sort | uniq > ${BED%.*}_geneids.txt
# Rscript 01_scripts/util/get_gene_names.R ${BED%.*}_geneids.txt

# echo "Summarizing results..."
# 01_scripts/util/get_genomic_context.R ${BED%.*}_dmr_context.txt
