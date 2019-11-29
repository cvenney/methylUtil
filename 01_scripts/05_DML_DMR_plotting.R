#!/usr/bin/env Rscript
## Methylation plotting

for (p in c("tidyverse", "data.table", "BiocManager", "ggplot2", "ggrepel", "ComplexHeatmap", "DSS", "bsseq", "vegan", "parallel", "GenomicAlignments")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("ComplexHeatmap", "DSS", "bsseq", "GenomicAlignments")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(library(p, character.only = T))}
    rm(p)
}

kyles_theme_adjustments <- theme_linedraw() + theme(axis.text = element_text(size = 12, colour = "black"), 
                                 axis.title = element_text(size = 14, colour = "black"),
                                 panel.grid = element_blank())

args <- commandArgs(T)
# args <- c("~/Projects/sasa_epi/sample_info_NC_027307.1.txt", "~/Projects/sasa_epi/03_results/DSS_two_group_test.txt.gz")

if (length(args) != 2)
    stop("Usage: DSS_DML_DMR.R <sample info file> <DM results>")

samples <- read.table(args[1], header = T, stringsAsFactors = FALSE)

if(!identical(colnames(samples), c("sample", "group", "family", "file")))
    stop("Samples file must contain a header row with names: \'sample\', \'group\', \'family\', \'file\'")


#### Load data ####

data_list <- lapply(1:nrow(samples), function(i) {
    file <- fread(samples[i, "file"], header = FALSE)[,c(-3:-4)]
    file <- file[V1 %like% c("NC_.*")]
    file[, V7 := V5 + V6]
    return(file[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)])
})

bs_obj <- makeBSseqData(data_list, samples[,"sample"])

rm(data_list)


#### Coverage filter ####

# Min and Max coverage

max_cov <- 20
min_cov <- 5

pass <- bs_obj@assays$data$Cov <= max_cov & bs_obj@assays$data$Cov >= min_cov
    
bs_obj <- bs_obj[rowSums(pass) >= 12,]    
rm(pass, max_cov, min_cov)

#### MDS plot of samples ####

un <- bs_obj@assays[[2]] - bs_obj@assays[[1]]
me <- bs_obj@assays[[2]]

M_values <- log2(me + 1) - log2(un + 1)
colnames(M_values) <- samples[, "sample"]
rm(un, me)

mds <- cmdscale(dist(t(M_values)), k = ncol(M_values) - 1)
colnames(mds) <- paste0("PC",1:ncol(mds))
mds <- data.frame(Sample = row.names(mds), Group = samples[, "group"], mds)

all_mds <- ggplot(mds, aes(x = PC1, y = PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_vline(xintercept = 0, colour = "grey") +
    kyles_theme_adjustments +
    geom_text_repel(aes(label = Sample), nudge_y = -10)
ggsave("03_results/global_methylation_mds_labs.png", plot = all_mds, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)

all_mds_points <- ggplot(mds, aes(x = PC1, y = PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_vline(xintercept = 0, colour = "grey") +
    geom_point(size = 2) +
    kyles_theme_adjustments
ggsave("03_results/global_methylation_mds_points.png", plot = all_mds_points, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)

rm(mds, all_mds, all_mds_points)


## Methylation ratios by sample

Beta_values <- bs_obj@assays$data$M / bs_obj@assays$data$Cov
# meth_ratio <- asin(sqrt(bs_obj@assays$data$M / bs_obj@assays$data$Cov))

df <- pivot_longer(as.data.frame(Beta_values), cols = starts_with("LJ"), names_to = "Sample", values_to = "Methylation")

methyl_ratio_hist <- ggplot(df, aes(x = Methylation)) +
    kyles_theme_adjustments +
    theme(axis.text = element_text(size = 8)) +
    geom_histogram(fill = "steelblue3", binwidth = 0.1) +
    facet_wrap(~Sample, ncol = 4) +
    ylab("Count")
ggsave("03_results/methylation_ratio_hist_by_sample.png", plot = methyl_ratio_hist, device = "png",
       width = 8, height = 8, units = "in", dpi = 300)
rm(df, methyl_ratio_hist, Beta_values)


#### RDA bewteen groups ####

data_rda <- t(as.matrix(M_values))
# data_rda <- t(as.matrix(Beta_values))

# Simple imputation of NaN with the mean methylation for each position
for (i in 1:ncol(data_rda)) {
    data_rda[is.nan(data_rda[,i]), i] <- mean(data_rda[,i], na.rm = TRUE)
}

parental_origin <- factor(rep(c("SAS", "Wild"), each = 8))

rda <- rda(data_rda ~ parental_origin)
fwrite(RsquareAdj(rda), "03_results/RDA_Rsquared.txt", sep = "\t")

# aov_results <- anova.cca(rda, parallel = 14)
# fwrite(aov_results, "03_results/RDA_ANOVA_results.txt", sep = "\t")

rm(data_rda, i)

sc <- scores(rda, choices = c(1:2), scaling = 3)
sc$sites <- data.frame(sc$sites, Pop = parental_origin)
rm(rda)

rda_sample_biplot <- ggplot(as.data.frame(sc$sites), aes(x = RDA1, y = PC1, color = Pop)) +
    kyles_theme_adjustments +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_point(size = 3) +
    #scale_color_manual(values = "") +
    geom_text(aes(label = rownames(as.data.frame(sc$sites))), nudge_x = 0.1, nudge_y = 0.1, size = 2)
ggsave("03_results/rda_sample_biplot.png", plot = rda_sample_biplot, device = "png",
       width = 4.5, height = 4, units = "in", dpi = 300)
rm(rda_sample_biplot)

rcov <- MASS::cov.rob(sc$species[, "RDA1"]) # robust covariance
resmaha <- mahalanobis(as.matrix(sc$species[, "RDA1"]), rcov$center, rcov$cov) # test statistic
lambda <- median(resmaha) / qchisq(0.5, df = 1)
sc$species <- as.data.frame(sc$species)
sc$species$pvals <- pchisq(resmaha / lambda, 1, lower.tail = FALSE)
rm(resmaha, lambda, rcov)
sc$species$fdr <- p.adjust(sc$species$pvals, method = "BH")
sc$species$chr <- as.factor(seqnames(bs_obj@rowRanges))
fwrite(sc$species, "03_results/rda_snp_scores.txt.gz", sep = "\t")

rda_CpG_biplot <- ggplot(sc$species, aes(x = RDA1, y = PC1, colour = fdr < 0.05)) +
    kyles_theme_adjustments +
    geom_point() +
    scale_color_manual(values = c("grey", "lightgreen"))
ggsave("03_results/rda_CpG_biplot.png", plot = rda_CpG_biplot, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)
rm(rda_CpG_biplot)

#rda_manhattan <- ggplot(sc$species, aes(x = 1:nrow(sc$species), y = abs(RDA1), colour = fdr < 0.05)) +
#    kyles_theme_adjustments +
#    geom_point(shape = "+")
#ggsave("03_results/rda_manhattan.png", plot = rda_manhattan, device = "png",
#       width = 8, height = 4, units = "in", dpi = 300)
rm(sc, rda_manhattan)


#### Diff Methyl results and heatmaps ####

dms <- fread(args[2])

all_density <- ggplot(dms, aes(x = diff * 100, fill = fdr < 0.05)) +
    geom_density() +
    kyles_theme_adjustments +
    scale_fill_brewer(type = "qual", palette = 4) +
    xlab("Change in methylation (%)") +
    ylab("Density")
ggsave("03_results/global_differential_methylation_density.png", plot = all_density, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)
rm(all_density)

sig_hist <- ggplot(dms[dms$fdr <0.05, ], aes(x = diff * 100)) +
    geom_histogram(binwidth = 2) +
    kyles_theme_adjustments +
    xlab("Change in methylation (%)") +
    ylab("Density")
ggsave("03_results/sig_differential_methylation_hist.png", plot = sig_hist, device = "png",
       width = 4.5, height = 4, units = "in", dpi = 300)
rm(sig_hist)


## Differential Methylation

dmls <- callDML(dms, delta = 0.2, p.threshold = 0.01)
fwrite(dmls, "03_results/differentially_methylated_loci_delta20_fdr0.01.txt", sep = "\t")
dmrs <- callDMR(dms, delta = 0.2, p.threshold = 0.01)
fwrite(dmrs, "03_results/differentially_methylated_regions_delta20_fdr0.01.txt", sep = "\t")

rm(dms, dmls)


## Heatmap

dmrs <- GRanges(
    seqnames = Rle(dmrs$chr),
    range = IRanges(start = dmrs$start, end = dmrs$end),
    nCG = dmrs$nCG,
    diffMethyl = dmrs$diff.Methy,
    areaStat = dmrs$areaStat
)
dmrs <- sort(dmrs)

Beta_values <- bs_obj@assays$data$M / bs_obj@assays$data$Cov
ME <- bs_obj@rowRanges
mcols(ME) <- as.data.frame(Beta_values)
# mcols(ME) <- as.data.frame(M_values)
rm(Beta_values)

hits <- findOverlaps(ME, dmrs, ignore.strand = TRUE)
dmr <- ME[queryHits(hits)]
mcols(dmrs) <- cbind(mcols(dmrs), aggregate(x = mcols(dmr), by = list(subjectHits(hits)), FUN = mean, na.rm = TRUE)[,-1])
rm(ME, dmr, hits)


png(filename = "03_results/all_DMR_heatmap.png", width = 8, height = 11, units = "in", res = 300)
Heatmap(
    matrix = as.matrix(mcols(dmrs)[samples$sample]),
    cluster_rows = TRUE,
    #row_split = factor(paste0(seqnames(dmr$dmrs), ":", start(dmr$dmrs), "-", end(dmr$dmrs))),
    row_title = NULL,
    #cluster_row_slices = TRUE,
    cluster_columns = FALSE,
    column_split = rep(c("SAS", "Wild"), each = 8),
    use_raster = TRUE,
    raster_device = "png"
)
dev.off()
