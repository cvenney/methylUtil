#!/usr/bin/env Rscript
## get_genomic_context.R

args <- commandArgs(TRUE)
# setwd("/Volumes/MacHD/BernatchezProjects/sasa_epi/methylUtil/"); args <- "06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001_dmr_context.txt"

if (length(args) != 1) {
    stop("Usage: get_genomic_context.R <bedtools_intersect_outfile>")
}

for (p in c("tidyverse", "data.table")) {
    if (!suppressPackageStartupMessages(require(p, character.only = T))) {
        install.packages(p, repos = "https://cloud.r-project.org", dependencies = T)
        suppressPackageStartupMessages(library(p, character.only = T))
    }
}


features <- c("region", "gene", "CDS", "exon", "mRNA", "three_prime_UTR", "five_prime_UTR", "transcript", "pseudogene", "lnc_RNA")

x <- fread(args[1], sep = "\t")

xtab <- as.data.frame(table(paste(x$V1, x$V2, sep = ":"), x$V9)) %>%
    pivot_wider(., id_cols = Var1, names_from = Var2, values_from = Freq)
xtab$feature_total <- rowSums(xtab[,-1])

if (any(!features %in% names(xtab))) {
    for (f in features[!features %in% names(xtab)]) {
        xtab <- mutate(xtab, !!f := 0)
    }
}

sign <- x %>% mutate(., Var1 = paste(V1, V2, sep = ":"), sign = sign(V5)) %>%
    select(., Var1, sign) %>% distinct()

context <- xtab %>%
    mutate(
        type = case_when(
            region == 1 & feature_total == 1 ~ "intergenic",
            gene > 0 & mRNA > 0 & exon == 0 ~ "intron",
            gene > 0 & transcript > 0 & exon == 0 ~ "intron",
            pseudogene > 0 ~ "pseudogene",
            gene > 0 & lnc_RNA > 0 ~ "lnc_RNA",
            gene > 0 & five_prime_UTR > 0 & three_prime_UTR == 0 ~ "five_prime_utr",
            gene > 0 & three_prime_UTR > 0 & five_prime_UTR == 0 ~ "three_prime_utr",
            gene > 0 & mRNA > 0 & exon > 0 ~ "exon",
            gene > 0 & transcript > 0 & exon > 0 ~ "exon",
            TRUE ~ "other"
            )
    )

context <- merge(context, sign, by = "Var1")

## Filter to investigate genomic contexts not captured by 
# context %>% filter(., type == "other")

output <- as.data.frame(table(context$type, context$sign)) %>%
    select("Feature" = Var1, "DM" = Var2, "N" = Freq) %>%
    pivot_wider(., names_from = DM, values_from = N)

write.table(output, paste0(sub("(.*)?\\..*$", "\\1", args[1]), "_summary.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

info <- unique.data.frame(x[,c(1:6)])
info$id <- paste(info$V1, info$V2, sep = ":")
context <- merge(info, context[, c("Var1", "type")], by.x = "id", by.y = "Var1")[,-1]
write.table(context, paste0(sub("(.*)?\\..*$", "\\1", args[1]), "_simple.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
