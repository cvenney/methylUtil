#!/usr/bin/env Rscript

args <- commandArgs(T)
# setwd("~/Projects/sasa_epi/methylUtil") ; args <- c("1000", "06_methylation_results/adults_6x6_min5_max20_groupWild_dmr_pval0.001.txt.gz", "06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001.txt.gz")

## Sanity checking
if (length(args) != 3)
    stop("Usage: DMR_enrichment.R <N reps> <DMR result> <BED file>")

## Load required packages
for (p in c("tidyverse", "parallel", "data.table", "BiocManager", "ggplot2", "cowplot", "GenomicRanges", "DSS")) {
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

if (any(grepl("^NC_02", idx$V1)))
    idx <- idx[V1 %like% "NC_02.*"]

max_bp <- idx[nrow(idx), V3]

# # Read and count DMC at pval and FDR thresholds
# cpg <- fread(sub("_dmr_.*", replacement = "_all_sites.txt.gz", args[2]))
# cpg <- cpg[!is.na(stat)]
# dmls <- cpg[pvals < as.numeric(sub(".*pval", "", sub(".txt.gz", "", args[2])))]
# fdrs <- cpg[fdrs < 0.05]
# cpg <- GRanges(
#     seqnames = cpg$chr,
#     ranges = IRanges(start = cpg$pos, width = 1),
#     strand = "+",
#     stat = cpg$stat,
#     pvals = 1
# )
# dmls <- GRanges(
#     seqnames = dmls$chr,
#     ranges = IRanges(start = dmls$pos, width = 1)
# )
# fdrs <- GRanges(
#     seqnames = fdrs$chr,
#     ranges = IRanges(start = fdrs$pos, width = 1)
# )

# Read and convert DMRs
dmrs <- fread(args[2])
dmrs <- GRanges(
    seqnames = Rle(dmrs$chr),
    ranges = IRanges(start = dmrs$start, end = dmrs$end)
)

# Read and convert gene features
features <- fread(args[3])
features <- GRanges(
    seqnames = Rle(features$chr),
    ranges = IRanges(start = features$start, end = features$end)
)

# calc real number of overlaps
# real_dml <- length(unique(queryHits(findOverlaps(dmls, features))))
# real_fdr <- length(unique(queryHits(findOverlaps(fdrs, features))))
real_dmr <- length(unique(queryHits(findOverlaps(dmrs, features))))

# Run resampling
sim <- mclapply(1:as.numeric(args[1]), mc.cores = detectCores(), function(s) {
    
    # resample regions based on sig. DMRs
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
    sim_dmrs1 <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = pos, width = width(dmrs))
    )
    chr1 <- sample(idx$V1, length(features), replace = TRUE, prob = idx$V3 / sum(idx$V3))
    pos1 <- numeric(length(features))
    for(i in 1:length(features)) {
        flag <- 1
        while(flag) {
            pos1[i] <- sample(idx[V1 == chr1[i], V2], 1)
            if ((pos1[i] + width(features[i])) <= idx[V1 == chr1[i], V2]) {
                flag <- 0
            }
        }
    }  
    sim_dmrs2 <- GRanges(
        seqnames = Rle(chr1),
        ranges = IRanges(start = pos1, width = width(features))
    )
    
    
    hits_dmr1 <- findOverlaps(sim_dmrs1, sim_dmrs2)
    sim_dmr1 <- length(unique(queryHits(hits_dmr1)))
    
    # # resample DMCs over all CpG sites
    # sim_sig_cpg <- sample(1:length(cpg), length(dmls), replace = FALSE)
    # sim_dmls <- cpg[sim_sig_cpg]
    # hits_dml <- findOverlaps(sim_dmls, features)
    # sim_dml <- length(unique(queryHits(hits_dml)))
    # 
    # # sub-sample for DMC passing FDR
    # sim_fdrs <- sim_dmls[sample(1:length(sim_dmls), length(fdrs), replace = FALSE)]
    # hits_fdr <- findOverlaps(sim_fdrs, features)
    # sim_fdr <- length(unique(queryHits(hits_fdr)))
    # 
    # # DMRs from random DMCs
    # dml_test <- cpg
    # dml_test[sim_sig_cpg]$pvals <- 0.0001
    # dml_test <- as.data.frame(dml_test)
    # colnames(dml_test)[1:2] <- c("chr", "pos")
    # class(dml_test) <- c(class(dml_test), "DMLtest.multiFactor")
    # dml_test <- suppressWarnings(callDMR(dml_test, delta = 0, p.threshold = 0.001))
    # if (is.null(dml_test)) {
    #     sim_dmr2 <- 0
    # } else {
    #     sim_dmrs2 <- GRanges(
    #         seqnames = Rle(dml_test$chr),
    #         ranges = IRanges(start = dml_test$start, end = dml_test$end)
    #     )
    #     hits_dmr2 <- findOverlaps(sim_dmrs2, features)
    #     sim_dmr2 <- length(unique(queryHits(hits_dmr2)))
    # }
    # 
    return(c(s, sim_dmr1))
})

df <- as.data.frame(do.call(rbind, sim))
colnames(df) <- c("Replicate", "DMR_region_overlaps")

# fwrite(df, sub("dmr_", "resampled_dmrs_", args[2]), sep = "\t", quote = FALSE)

# dml_hist <- ggplot(df, aes(x = DMC_overlaps)) +
#     geom_histogram(binwidth = 5) +
#     theme_classic() +
#     labs(x = "Number of Overlapped Genes", y = "Frequency", title = paste0("DMC: N = ", real_dml, "; p = ", sum(df$DMC_overlaps > real_dml) / nrow(df))) +
#     theme(axis.text = element_text(size = 8, colour = "black"),
#           axis.title = element_text(size = 8, colour = "black"),
#           plot.title = element_text(size = 10, face = "bold")) +
#     geom_vline(xintercept = real_dml, colour = "red")
# 
# fdr_hist <- ggplot(df, aes(x = FDR_overlaps)) +
#     geom_histogram(binwidth = 5) +
#     theme_classic() +
#     labs(x = "Number of Overlapped Genes", y = "Frequency", title = paste0("DMC FDR: N = ", real_fdr, "; p = ", sum(df$FDR_overlaps > real_fdr) / nrow(df))) +
#     theme(axis.text = element_text(size = 8, colour = "black"),
#           axis.title = element_text(size = 8, colour = "black"),
#           plot.title = element_text(size = 10, face = "bold")) +
#     geom_vline(xintercept = real_fdr, colour = "red")

dmr1_hist <- ggplot(df, aes(x = DMR_region_overlaps)) +
    geom_histogram(binwidth = 5) +
    theme_classic() +
    labs(x = "Number of Overlapped Genes", y = "Frequency", title = paste0("DMR region: N = ", real_dmr, "; p = ", (df$DMR_region_overlaps > real_dmr) / nrow(df))) +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 8, colour = "black"),
          plot.title = element_text(size = 10, face = "bold")) +
    geom_vline(xintercept = real_dmr, colour = "red")

dmr1_hist

# dmr2_hist <- ggplot(df, aes(x = DMR_from_DMC_overlaps)) +
#     geom_histogram(binwidth = 5) +
#     theme_classic() +
#     labs(x = "Number of Overlapped Genes", y = "Frequency", title = paste0("DMR <- DMC: N = ", real_dmr, "; p = ", (df$DMR_from_DMC_overlaps > real_dmr) / nrow(df))) +
#     theme(axis.text = element_text(size = 8, colour = "black"),
#           axis.title = element_text(size = 8, colour = "black"),
#           plot.title = element_text(size = 10, face = "bold")) +
#     geom_vline(xintercept = real_dmr, colour = "red")

# save_plot(sub(".txt.gz", "_resampled_gene_overlap.png", args[2]), plot_grid(dml_hist, fdr_hist, dmr1_hist, dmr2_hist, ncol = 2, nrow = 2), base_height = 7, base_width = 7)
