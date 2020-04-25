#!/usr/bin/env Rscript

args <- commandArgs(T)
# setwd("/Volumes/MacHD/BernatchezProjects/sasa_epi/") ; args <- c("1000", "06_methylation_results/adults_6x7_min5_max20_groupWild_dmr_pval0.001.txt.gz", "05_bed_files/gene_w_promotors.bed")

## Sanity checking
if (length(args) != 3)
    stop("Usage: DMR_enrichment.R <N reps> <DMR result> <BED file>")

## Load required packages
for (p in c("tidyverse", "data.table", "BiocManager", "ggplot2", "GenomicRanges")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("ComplexHeatmap", "DSS", "bsseq", "GenomicAlignments")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(library(p, character.only = T))}
    rm(p)
}

if (!file.exists("02_reference/genome.fasta.fai"))
    stop("Genome index file 'genome.fasta.fai' not found!")

idx <- fread("02_reference/genome.fasta.fai")
idx <- idx[V1 %like% "NC_*"]

max_bp <- idx[nrow(idx), V3]

dmrs <- fread(args[2])
dmrs <- GRanges(
    seqnames = Rle(dmrs$chr),
    ranges = IRanges(start = dmrs$start, end = dmrs$end)
)

features <- fread(args[3])
features <- GRanges(
    seqnames = Rle(features$V1),
    ranges = IRanges(start = features$V2, end = features$V3)
)

real <- length(unique(queryHits(findOverlaps(dmrs, features))))

sim <- numeric(length = as.numeric(args[1]))

prog_threshold <- as.numeric(args[1])/(as.numeric(args[1])*0.05)

for (s in 1:length(sim)) {
    
    if( s %% prog_threshold == 0)
        message(paste0("Completed:", s/prog_threshold, "% of samplings"))
    chr <- sample(idx$V1, length(dmrs), replace = TRUE, prob = idx$V3 / sum(idx$V3))
    pos <- numeric(length(dmrs))
    for(i in 1:length(dmrs)) {
        flag <- 1
        while(flag) {
            pos[i] <- sample(idx[V1 == chr[i], V2], 1)
            if ((pos[i] + width(dmrs[i])) <= idx[V1 == chr[i], V2]) {
                flag <- 0
            }
        }
    }  
    sim_dmrs <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = pos, width = width(dmrs))
    )
    
    hits <- findOverlaps(sim_dmrs, features)
    
    sim[s] <- length(unique(queryHits(hits)))
}

hist(sim)
abline(v = real)

sum(sim >= real) / length(sim)
