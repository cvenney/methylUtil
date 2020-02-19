#!/bin/bash
#get_geneids_for_dmr.sh

if [[ $# -ne 2 ]]; then
	echo "USAGE: $0 <dmr_context> <gene/trans>"
	exit
fi

DMR=$1
TYPE=$2

if [[ $TYPE == "gene" ]] 
then
	awk '($9 == "gene") {print $15}' $DMR | \
	sed -e 's/ID=//' -e 's/;.*GeneID:/:/' -e 's/;.*//' |
	sort | uniq
elif [[ $TYPE == "trans" ]]
then
	echo "Currently not implemented"
	exit
else
	echo "ERROR: option \'$TYPE\' not supported"
	exit
fi

