#!/bin/bash
#convert_DMR_txt_to_bed.sh

SUBSETS="07_jackknife_results"

for i in $(ls -1 "$SUBSETS"/*.txt | perl -pe 's/\.txt//g')
do
	echo $(basename $i)
	tail -n +2 "$i".txt | tr -d \ |
#	awk 'NR == 1 {print "chrom", "start", "end"} NR > 1 {gsub(/"/, "", $2); print $2, $3, $4}' OFS='\t' "$i".txt > "$i".bed
	awk 'NR > 1 {gsub(/"/, "", $2); print $2, $3, $4}' OFS='\t' "$i".txt > "$i".bed
done