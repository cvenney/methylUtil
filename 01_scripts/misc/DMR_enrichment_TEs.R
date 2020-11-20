#!/usr/bin/env Rscript

args <- commandArgs(T)
# setwd("~/Projects/sasa_epi/methylUtil") ; args <- c("1000", "06_methylation_results/adults_6x6_min5_max20_groupWild_dmr_pval0.001.txt.gz", "05_bed_files/genes.bed")

setwd("~/Projects/sasa_epi/methylUtil") ; args <- c("1000", "06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001.txt.gz", "05_bed_files/TEs.bed")

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

# Read and count DMC at pval and FDR thresholds
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
    seqnames = Rle(features$V1),
    ranges = IRanges(start = features$V2, end = features$V3),
    type = features$V4
)

del <- features[!features$type %like% "Unknown"]
# dup <- features[features$type == "<DUP>"]
# inv <- features[features$type == "<INV>"]

# calc real number of overlaps
real_del <- length(unique(subjectHits(findOverlaps(dmrs, del))))
# real_dup <- length(unique(subjectHits(findOverlaps(dmrs, dup))))
# real_inv <- length(unique(subjectHits(findOverlaps(dmrs, inv))))

# Run resampling
sim <- mclapply(1:as.numeric(args[1]), mc.cores = detectCores(), function(s) {
    
    # resample regions based on sig. DMRs
    chr <- sample(idx$V1, length(dmrs), replace = TRUE, prob = idx$V2 / sum(idx$V2))
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
    hits_del <- length(unique(subjectHits(findOverlaps(sim_dmrs, del))))
    # hits_dup <- length(unique(subjectHits(findOverlaps(sim_dmrs, dup))))
    # hits_inv <- length(unique(subjectHits(findOverlaps(sim_dmrs, inv))))
    
    return(c(s, hits_del))
})

df <- as.data.frame(do.call(rbind, sim))
# colnames(df) <- c("Replicate", "DMC_overlaps", "FDR_overlaps", "DMR_region_overlaps", "DMR_from_DMC_overlaps")
colnames(df) <- c("Replicate", "TE_overlaps")

fwrite(df, sub("dmr_", "dmr_TE_overlaps", args[2]), sep = "\t", quote = FALSE)

del_hist <- ggplot(df, aes(x = TE_overlaps)) +
    geom_histogram(binwidth = 1) +
    theme_classic() +
    labs(x = "Number of Overlapped TEs", y = "Frequency", title = paste0("True Overlaps: N = ", real_del, "; p = ", sprintf("%.3f", as.numeric(sum(df$TE_overlaps > real_del) / nrow(df))))) +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 8, colour = "black"),
          plot.title = element_text(size = 10, face = "bold")) +
    geom_vline(xintercept = real_del, colour = "red")

# dup_hist <- ggplot(df, aes(x = duplication_overlaps)) +
#     geom_histogram(binwidth = 1) +
#     theme_classic() +
#     labs(x = "Number of Overlapped Duplications", y = "Frequency", title = paste0("True Overlaps: N = ", real_dup, "; p = ", sum(df$duplication_overlaps > real_dup) / nrow(df))) +
#     theme(axis.text = element_text(size = 8, colour = "black"),
#           axis.title = element_text(size = 8, colour = "black"),
#           plot.title = element_text(size = 10, face = "bold")) +
#     geom_vline(xintercept = real_dup, colour = "red")
# 
# inv_hist <- ggplot(df, aes(x = inversion_overlaps)) +
#     geom_histogram(binwidth = 1) +
#     theme_classic() +
#     labs(x = "Number of Overlapped Inversions", y = "Frequency", title = paste0("True Overlaps: N = ", real_inv, "; p = ", sum(df$inversion_overlaps > real_inv) / nrow(df))) +
#     theme(axis.text = element_text(size = 8, colour = "black"),
#           axis.title = element_text(size = 8, colour = "black"),
#           plot.title = element_text(size = 10, face = "bold")) +
#     geom_vline(xintercept = real_inv, colour = "red")


save_plot(sub(".txt.gz", "_resampled_SV_overlap.png", args[2]), del_hist, base_height = 3.5, base_width = 8.5)
