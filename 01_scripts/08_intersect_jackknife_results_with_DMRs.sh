#!/bin/bash
#intersect_jackknife_results_with_DMRs.sh

FULL_AD_DMRs="06_methylation_results/DSS_results_min5_max80_adult_tempfroid_dmr_pval0.001_noheader.txt"
FULL_JU_DMRs="06_methylation_results/DSS_results_min5_max80_juv_tempfroid_dmr_pval0.001_noheader.txt"
FULL_AD_JU_DMRs="06_methylation_results/DSS_results_min5_max80_adult_tempfroid.juv_tempfroid_dmr_pval0.001_noheader.txt"
SUBSETS="07_jackknife_results"

module load bedtools/2.30.0	

# adult
bedtools intersect -a "$FULL_AD_DMRs" -b "$SUBSETS"/brook_charr_all_CpGs_F*adult_tempfroid_DMRs*.bed -C -filenames -f 0.8 > "$SUBSETS"/DMR_jackknife_overlap_counts_80pctoverlap_adulttemp.txt

# juv
bedtools intersect -a "$FULL_JU_DMRs" -b "$SUBSETS"/brook_charr_all_CpGs_F*C_juv_tempfroid_DMRs*.bed "$SUBSETS"/brook_charr_all_CpGs_F*F_juv_tempfroid_DMRs*.bed -C -filenames -f 0.8 > "$SUBSETS"/DMR_jackknife_overlap_counts_80pctoverlap_juvtemp.txt

# adult x juv
bedtools intersect -a "$FULL_AD_JU_DMRs" -b "$SUBSETS"/brook_charr_all_CpGs_F*adult_tempfroid_juv_tempfroid_DMRs*.bed -C -filenames -f 0.8 > "$SUBSETS"/DMR_jackknife_overlap_counts_80pctoverlap_adultxjuvtemp.txt