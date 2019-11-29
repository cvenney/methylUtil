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

kyles_theme_adjustments <- theme_linedraw() + theme(axis.text = element_text(size = 12, colour = "black"), 
                                                    axis.title = element_text(size = 14, colour = "black"),
                                                    panel.grid = element_blank())

args <- commandArgs(T)
# args <- c(1, "~/Projects/safo_epi/02_data/sample_info.txt", "~/Projects/safo_epi/02_data/chrs.txt")

## Sanity checking
if (length(args) != 3)
    stop("Usage: DSS_coverage_plotting.R <n_core> <sample_info.txt> <chrs.txt>")

n_cores <- as.integer(args[1])
setDTthreads(threads = n_cores)
samples <- read.table(args[2], header = T, stringsAsFactors = FALSE)
chrs <- read.table(args[3], header = FALSE, stringsAsFactors = FALSE)[,1]

ncol <- floor(sqrt(length(chrs)))

if(!all(c("sample", "file") %in% colnames(samples)))
    stop("Samples file must contain a header row with names: \'sample\', \'file\'.")


## Read in data chr by chr
bs_obj_all <- lapply(chrs, function(chr) {
    
    # Create local sample info and adujust filenames for specific chr
    lsamples <- samples
    lsamples$file <- sub("\\.bedGraph\\.gz", paste0("_", chr, "\\.bedGraph\\.gz"), lsamples$file)
    
    message(paste0("Loading chr: ", chr))
    
    # Load data and convert to BSseq object
    data_list <- lapply(1:nrow(samples), function(i) {
        file <- fread(lsamples[i, "file"], header = FALSE)[,c(-3:-4)]
        file[, V7 := V5 + V6]
        return(file[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)])
    })
    bs_obj <- makeBSseqData(data_list, samples[,"sample"])
    rm(data_list)
    return(bs_obj)
})

bs_obj_all <- suppressWarnings(do.call(rbind, bs_obj_all))

chr <- as.character(seqnames(bs_obj_all))

df_cpg <- data.frame(chr = as.factor(if_else(chr %in% chrs, chr, "contigs")),
                     median = rowMedians(getCoverage(bs_obj_all, type = "Cov")),
                     mean = rowMeans(getCoverage(bs_obj_all, type = "Cov")))
rm(chr)

median_hist <- ggplot(df_cpg, aes(x = median)) +
    theme_linedraw() + 
    theme(panel.grid = element_blank(), 
          axis.text.y = element_blank(),
          legend.position = c(0.9, 0.1)) +
    geom_histogram(binwidth = 2) +
    xlim(c(-0.1,40)) +
    geom_vline(data = data.frame(quantile = factor(c("5X", "10X", "0.99", "0.995"), levels = c("5X", "10X", "0.99", "0.995")), 
                                 value = c(5, 10, quantile(df_cpg$mean, c(0.99, 0.995)))),
               aes(xintercept = value, linetype = quantile)) +
    facet_wrap(~ chr, ncol = ncol, scales = "free_y")

suppressWarnings(ggsave("06_methylation_results/coverage/median_coverage_hist.png", median_hist, device = "png", width = 8.5, height = 11, units = "in", dpi = 300))
rm(median_hist)

mean_hist <- ggplot(df_cpg, aes(x = mean)) +
    theme_linedraw() + 
    theme(panel.grid = element_blank(), 
          axis.text.y = element_blank(),
          legend.position = c(0.9, 0.1)) +
    geom_histogram(binwidth = 2) +
    xlim(c(-0.1,40)) +
    geom_vline(data = data.frame(quantile = factor(c("5X", "10X", "0.99", "0.995"), levels = c("5X", "10X", "0.99", "0.995")), 
                                 value = c(5, 10, quantile(df_cpg$mean, c(0.99, 0.995)))),
               aes(xintercept = value, linetype = quantile)) +
    facet_wrap(~ chr, ncol = ncol, scales = "free_y")

suppressWarnings(ggsave("06_methylation_results/coverage/mean_coverage_hist.png", mean_hist, device = "png", width = 8.5, height = 11, units = "in", dpi = 300))
rm(mean_hist, df_cpg)


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
write.table(ind_coverage, "06_methylation_results/coverage/Individual_coverage_stats.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


## Principle coordinates analysis to identify outlier individuals based on coverage
euc_dist <- dist(t(getCoverage(bs_obj_all, type = "Cov")))
mds <- as.data.frame(cmdscale(euc_dist))
names(mds) <- c("PC1", "PC2")
mds$Ind <- row.names(mds)
mds_plot <- ggplot(mds, aes(x = PC1, y = PC2, label = Ind)) +
    kyles_theme_adjustments +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_label_repel(fill = "black", color = "white")
ggsave("06_methylation_results/coverage/mds_plot_coverage.png", mds_plot, device = "png", width = 8, height = 8, units = "in", dpi = 300)

euc_miss <- dist(t(getCoverage(bs_obj_all, type = "Cov")>0))
mds <- as.data.frame(cmdscale(euc_miss))
names(mds) <- c("PC1", "PC2")
mds$Ind <- row.names(mds)
mds_plot <- ggplot(mds, aes(x = PC1, y = PC2, label = Ind)) +
    kyles_theme_adjustments +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_label_repel(fill = "black", color = "white")
ggsave("06_methylation_results/coverage/mds_plot_missing_data.png", mds_plot, device = "png", width = 8, height = 8, units = "in", dpi = 300)
