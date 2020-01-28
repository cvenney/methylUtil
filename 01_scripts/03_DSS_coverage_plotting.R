#!/usr/bin/env Rscript
## Script for generating plot to inform coverage thresholds
## for DML / DMR quantification using DSS

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "DSS", "bsseq", "parallel", "tidyverse", "ggrepel")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("DSS", "bsseq")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

theme_adjustments <- theme_linedraw() + theme(axis.text = element_text(size = 12, colour = "black"), 
                                              axis.title = element_text(size = 14, colour = "black"),
                                              panel.grid = element_blank())

args <- commandArgs(T)
# args <- c(8, "~/Projects/safo_epi/methylUtil/sample_info_unpaired_control.txt"); setwd("~/Projects/safo_epi/methylUtil")
# args <- c(8, "~/Projects/sasa_epi/methylUtil/juvenile_samples_8x8.txt"); setwd("~/Projects/sasa_epi/methylUtil")

## Sanity checking
if (length(args) != 3)
    stop("Usage: DSS_coverage_plotting.R <n_core> <sample_info.txt> <chrs.txt>")

n_cores <- as.integer(args[1])
setDTthreads(threads = n_cores)
samples <- read.table(args[2], header = T, stringsAsFactors = FALSE)

ncol <- floor(sqrt(length(chrs)))

if(!all(c("sample", "file") %in% colnames(samples)))
    stop("Samples file must contain a header row with names: \'sample\', \'file\'.")

outfile_name <- paste0("06_methylation_results/coverage/", gsub("\\..*", "", basename(args[2])))

## Read in data
bs_obj_all <- lapply(1:nrow(samples), function(i) {
    message(paste0("Loading sample: ", samples[i, "sample"]))
    samp <- fread(samples[i, "file"], header = FALSE)[,c(-3:-4)]
    samp[, V7 := V5 + V6] # Combine me+ and me- counts for total coverage
    #bs_obj <- makeBSseqData(list(samp[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)]), samples[i,"sample"])
    return(samp[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)])
})

bs_obj_all <- makeBSseqData(dat = bs_obj_all, sampleNames = samples$sample)
bs_obj_all <- bs_obj_all[(rowSums(getCoverage(bs_obj_all, type = "Cov") >= 1) == ncol(bs_obj_all)), ]

saveRDS(bs_obj_all, paste0("06_methylation_results/", gsub("\\..*", "", basename(args[2])), "_all_data.rds"), compress = "gzip")


## Generate summary table of coverage statistics by individual
ind_coverage <- data.frame(
    Ind = samples$sample,
    Total_Coverage = colSums(getCoverage(bs_obj_all, type = "Cov")),
    N_CpG_lt_eq_0 = colSums(getCoverage(bs_obj_all, type = "Cov")>=0),
    N_CpG_lt_eq_5 = colSums(getCoverage(bs_obj_all, type = "Cov")>=5),
    N_CpG_lt_eq_10 = colSums(getCoverage(bs_obj_all, type = "Cov")>=10),
    Median_Coverage = colMedians(getCoverage(bs_obj_all, type = "Cov")),
    Mean_Coverage = colMeans(getCoverage(bs_obj_all, type = "Cov")),
    SD_Coverage = colSds(getCoverage(bs_obj_all, type = "Cov")),
    Max_Coverage = colMaxs(getCoverage(bs_obj_all, type = "Cov"))
)
write.table(ind_coverage, paste0(outfile_name, "_Individual_coverage_stats.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


## Principle coordinates analysis to identify outlier individuals based on coverage
euc_dist <- dist(t(getCoverage(bs_obj_all, type = "Cov")))
mds <- as.data.frame(cmdscale(euc_dist))
names(mds) <- c("PC1", "PC2")
mds$Ind <- row.names(mds)
mds_plot <- ggplot(mds, aes(x = PC1, y = PC2, label = Ind)) +
    theme_adjustments +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_label_repel(fill = "black", color = "white")
ggsave(paste0(outfile_name, "_mds_plot_coverage.png"), mds_plot, device = "png", width = 8, height = 8, units = "in", dpi = 300)

euc_miss <- dist(t(getCoverage(bs_obj_all, type = "Cov")>0))
mds <- as.data.frame(cmdscale(euc_miss))
names(mds) <- c("PC1", "PC2")
mds$Ind <- row.names(mds)
mds_plot <- ggplot(mds, aes(x = PC1, y = PC2, label = Ind)) +
    theme_adjustments +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_label_repel(fill = "black", color = "white")
ggsave(paste0(outfile_name, "_mds_plot_missing_data.png"), mds_plot, device = "png", width = 8, height = 8, units = "in", dpi = 300)
