#!/usr/bin/env Rscript
## get_genomic_context.R

if (suppressPackageStartupMessages(!require(tidyverse))) {
    install.packages("tidyverse");suppressPackageStartupMessages(require(tidyverse))
}

args <- commandArgs(TRUE)
# setwd("~/Desktop/salmon_annotation/"); args <- "dmr_context.txt"

if (length(args) != 1) {
    stop("Usage: get_genomic_context.R <bedtools_intersect_outfile>")
}

x <- read.table(args[1], sep = "\t", stringsAsFactors = FALSE)

xtab <- as.data.frame(table(paste(x$V1, x$V2, sep = ":"), x$V9)) %>%
    pivot_wider(., id_cols = Var1, names_from = Var2, values_from = Freq)
xtab$feature_total <- rowSums(xtab[,-1])

context <- xtab %>%
    mutate(
        type = case_when(
            region == 1 & feature_total == 1 ~ "intergenic",
            gene > 0 & mRNA > 0 & exon == 0 ~ "intron",
            pseudogene > 0 ~ "pseudogene",
            gene > 0 & lnc_RNA > 0 ~ "lnc_RNA",
            gene > 0 & five_prime_UTR > 0 & three_prime_UTR == 0 ~ "five_prime_utr",
            gene > 0 & three_prime_UTR > 0 & five_prime_UTR == 0 ~ "three_prime_utr",
            gene > 0 & mRNA > 0 & exon > 0 ~ "exon",
            TRUE ~ "other"
            )
    )

## Filter to investigate genomic contexts not captured by 
# context %>% filter(., type == "other")

output <- as.data.frame(table(context$type)) %>%
    select("Feature" = Var1, "N" = Freq)

write.table(output, paste0("06_methylation_results/", sub("05_bed_files/", "", args[1])),
            row.names = FALSE, quote = FALSE, sep = "\t")