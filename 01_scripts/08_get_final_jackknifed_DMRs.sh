#!/bin/bash
#get_final_jackknifed_DMRs.sh

AD_OVERLAPS="07_jackknife_results/DMR_jackknife_overlap_counts_80pctoverlap_adulttemp.txt"
JU_OVERLAPS="07_jackknife_results/DMR_jackknife_overlap_counts_80pctoverlap_juvtemp.txt"
AD_JU_OVERLAPS="07_jackknife_results/DMR_jackknife_overlap_counts_80pctoverlap_adultxjuvtemp.txt"

### adult
#   require 8 overlaps (from 8 unique subsets) for DMR to be considered verified
awk -F"\t" '$8>=1' "$AD_OVERLAPS" |
#cut -f -3 | uniq -c | awk '$1>=8' > 07_jackknife_results/DMR_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min8.txt  # for maternal jackknifing
cut -f -3 | uniq -c | awk '$1>=14' > 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adulttemp.txt  # for family jackknifing

awk 'NR {print $2, $3, $4}' OFS='\t' \
	07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adulttemp.txt \
	> 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adulttemp.bed



### juvenile
#   require 8 overlaps (from 8 unique subsets) for DMR to be considered verified
awk -F"\t" '$8>=1' "$JU_OVERLAPS" |
#cut -f -3 | uniq -c | awk '$1>=8' > 07_jackknife_results/DMR_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min8.txt  # for maternal jackknifing
cut -f -3 | uniq -c | awk '$1>=14' > 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_juvtemp.txt  # for family jackknifing

awk 'NR {print $2, $3, $4}' OFS='\t' \
	07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_juvtemp.txt \
	> 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_juvtemp.bed



### adult x juvenile
#   require 8 overlaps (from 8 unique subsets) for DMR to be considered verified
awk -F"\t" '$8>=1' "$AD_JU_OVERLAPS" |
#cut -f -3 | uniq -c | awk '$1>=8' > 07_jackknife_results/DMR_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min8.txt  # for maternal jackknifing
cut -f -3 | uniq -c | awk '$1>=14' > 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adultxjuvtemp.txt  # for family jackknifing

awk 'NR {print $2, $3, $4}' OFS='\t' \
	07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adultxjuvtemp.txt \
	> 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adultxjuvtemp.bed

# get beta values for heatmap
#module load bedtools
# need to hash out first row of -a
#bedtools intersect -header -f 0.8 \
#	-a 06_methylation_results_main/DSS_results_adult_temp_DMR_heatmap_mean_Betas.txt \
#	-b 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adulttemp.bed \
#	> 07_jackknife_results/DMR_family_jackknife_overlap_counts_80pctoverlap_only_covered_with_count_min14_adulttemp_heatmap_mean_Betas_f0.8.txt