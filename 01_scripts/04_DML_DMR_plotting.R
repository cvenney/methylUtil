#!/usr/bin/env Rscript
## Methylation plotting

for (p in c("configr", "png", "tidyverse", "data.table", "BiocManager", "ggplot2", "ggrepel", "ComplexHeatmap", "DSS", "bsseq", "vegan", "parallel")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("ComplexHeatmap", "DSS", "bsseq", "GenomicAlignments")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(library(p, character.only = T))}
    rm(p)
}

theme_adjustments <- theme_linedraw() + theme(axis.text = element_text(size = 12, colour = "black"), 
                                              axis.title = element_text(size = 14, colour = "black"),
                                              panel.grid = element_blank())

args <- commandArgs(T)
# args <- "~/Projects/sasa_epi/methylUtil/config_juveniles_7x7.yml" ; setwd("~/Projects/sasa_epi/methylUtil")
# args <- "config_unpaired.yml" ; setwd("~/Projects/safo_epi/methylUtil")

## Sanity checking
if (length(args) != 1)
    stop("Usage: 05_DML_DMR_plotting.R <config.yml>")

if (!is.yaml.file(args[1]))
    stop("You must supply a configuration file in YAML format.\nUsage: 05_DML_DMR_plotting.R <config.yml>")

config <- read.config(args[1])


# Parse formula
if (is.null(config$options$formula) | !grepl("\\~", config$options$formula))
    stop("Invalid formula. You must provide a design formula beginning with a tilde (e.g. \'~ Treatment\')")

formula <- as.formula(config$options$formula)
formula_parts <- unlist(strsplit(config$options$formula, split = c("\\~ |\\~|\\s\\+\\s|\\s\\+|\\+\\s|\\+|\\*|\\s\\*|\\s\\*\\s|\\*\\s|\\:")))[-1]


# Load files and set up design
samples <- read.table(config$input$sample_info, header = T, stringsAsFactors = FALSE)

if(!all(c("sample", "file", formula_parts) %in% colnames(samples)))
    stop("Samples file must contain a header row with names: \'sample\', \'file\', and given factor(s).")

if (any(!formula_parts %in% colnames(samples)))
    stop("Factors specified in formula design are not present in the ")
design <- data.frame(samples[, formula_parts])
design[] <- lapply(design, factor)
names(design) <- formula_parts


# Set number of cores for parallel computing
if (is.null(config$options$n_cores)) {
    warning("\'n_cores\' not specified. Default to using 1 core.")
    n_cores <- 1
} else {
    if (config$options$n_cores == 0)
        warning("Using all cores will require a lot of memory")
    n_cores <- ifelse(config$options$n_cores == 0, detectCores(), config$options$n_cores)
    setDTthreads(n_cores)
}

# Set coverage options for filtering
if (is.null(config$options$max_coverage)) {
    warning("\'max_coverage\' not specified. Default to using a value of 30.")
    max_cov <- 30L
} else {
    max_cov <- config$options$max_coverage
}

if (is.null(config$options$max_coverage)) {
    warning("\'max_coverage\' not specified. Default to using a value of 10.")
    min_cov <- 10L
} else {
    min_cov <- config$options$min_coverage
}

if (is.null(config$options$max_coverage)) {
    warning("\'min_individuals\' not specified. Requiring coverage for all individuals.")
    min_ind <- nrow(samples)
} else {
    min_ind <- config$options$min_individuals
}

if (is.null(config$options$pval_threshold) | !is.numeric(config$options$pval_threshold)) {
    warning("Invalid pval threshold. Default to using a value of 1e-5.")
    pval <- 1e-5
} else {
    pval <- as.numeric(config$options$pval_threshold)
}

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE) & (is.null(config$options$delta) | !is.numeric(config$options$delta))) {
    warning("Delta required for Wald tests. Default to using a value of 0.1.")
    delta <- 0.1
} else {
    delta <- config$options$delta
}

bs_obj_path <- paste0(config$output$outfile_prefix, "_min", min_cov, "_max", max_cov)

if (file.exists(bs_obj_path)) {
    bs_obj_all <- readRDS(file = bs_obj_path)
} else {
    stop("BSseq object not found! You need to run 03_DSS_model.R first.")
}


#### MDS plot of samples ####

un <- getCoverage(bs_obj_all, type = "Cov") - getCoverage(bs_obj_all, type = "M")
me <- getCoverage(bs_obj_all, type = "M")

M_values <- log2(me + 1) - log2(un + 1)
colnames(M_values) <- samples[, "sample"]
rm(un, me)

mds <- cmdscale(dist(t(M_values)), k = ncol(M_values) - 1)
colnames(mds) <- paste0("PC",1:ncol(mds))
mds <- data.frame(Sample = row.names(mds), Group = interaction(samples[,formula_parts]), mds)
fwrite(mds, paste0(config$output$outfile_prefix, "_methylation_MDS.txt"), quote = FALSE, sep = "\t")

all_mds <- ggplot(mds, aes(x = PC1, y = PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_vline(xintercept = 0, colour = "grey") +
    theme_adjustments +
    geom_text_repel(aes(label = Sample), nudge_y = -10)
ggsave(paste0(config$output$outfile_prefix, "_global_methylation_mds_labs.png"), plot = all_mds, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)

all_mds_points <- ggplot(mds, aes(x = PC1, y = PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_vline(xintercept = 0, colour = "grey") +
    geom_point(size = 2) +
    theme_adjustments
ggsave(paste0(config$output$outfile_prefix, "_global_methylation_mds_points.png"), plot = all_mds_points, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)

rm(mds, all_mds, all_mds_points, M_values)


## Methylation ratios by sample

Beta_values <- getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")
# Beta_values <- asin(sqrt(getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")))

df <- pivot_longer(as.data.frame(Beta_values), everything(), names_to = "Sample", values_to = "Methylation")
rm(Beta_values)
Beta_summary <- df %>% 
    group_by(Sample) %>% 
    summarise(Mean = mean(Methylation, na.rm = TRUE), 
              Median = median(Methylation, na.rm = TRUE), 
              SD = sd(Methylation, na.rm = TRUE), 
              NAs = sum(is.na(Methylation)))

fwrite(Beta_summary, paste0(config$output$outfile_prefix, "_Beta_summary_by_individual.txt"), quote = FALSE, sep = "\t")

methyl_ratio_hist <- ggplot(df, aes(x = Methylation)) +
    theme_adjustments +
    theme(axis.text = element_text(size = 8)) +
    geom_histogram(fill = "steelblue3", binwidth = 0.05) +
    facet_wrap(~Sample, ncol = floor(sqrt(nrow(samples)))) +
    ylab("Count")
ggsave(paste0(config$output$outfile_prefix, "_methylation_ratio_hist_by_sample.png"), plot = methyl_ratio_hist, device = "png",
       width = 8.5, height = 11, units = "in", dpi = 300)
rm(df, methyl_ratio_hist)



#### Diff Methyl results and heatmaps ####
dml_dmr_summary <- function(dmls, dmrs, coef = NULL, flag = NULL) {
    
    if (flag == 0 ) {
        dml_summary <- dmls[, .("Diff_Mean" = mean(diff), 
                                "Diff_Median" = median(diff),
                                "Diff_SD" = sd(diff),
                                "stat_Mean" = mean(stat), 
                                "stat_Median" = median(stat),
                                "stat_SD" = sd(stat),
                                "N" = .N), by = sign(stat)]
    } else if (flag == 1) {
        dml_summary <- dmls[, .("stat_Mean" = mean(stat), 
                                "stat_Median" = median(stat),
                                "stat_SD" = sd(stat),
                                "N" = .N), by = sign(stat)]
    }
    
    fwrite(dml_summary, paste0(config$output$outfile_prefix, "_DML_hyper_hypo_distribution_", coef, ".txt"), quote = FALSE, sep = "\t")
    
    dmr_summary <- dmrs[, .("areaStat_Mean" = mean(areaStat), 
                            "areaStat_Median" = median(areaStat),
                            "areaStat_SD" = sd(areaStat),
                            "nCG" = sum(nCG),
                            "nCG_mean" = mean(nCG),
                            "nCG_median" = median(nCG),
                            "nCG_SD" = sd(nCG),
                            "N" = .N), by = sign(areaStat)]

    # "length_Mean" = mean(end - start),
    # "length_Median" = median(end - start),
    # "length_SD" = sd(end - start),
    
    fwrite(dmr_summary, paste0(config$output$outfile_prefix, "_DMR_hyper_hypo_distribution_", coef, ".txt"), quote = FALSE, sep = "\t")

}

DMR_heatmap <- function(dmrs, Betas, design, anno_columns, sample_info, coef = NULL) {
    
    if(!"GRanges" %in% class(dmrs)) {
        ldmrs <- GRanges(
            seqnames = Rle(dmrs$chr),
            range = IRanges(start = dmrs$start, end = dmrs$end)
        )
        ldmrs <- sort(ldmrs)
    } else {
        ldmrs <- dmrs
    }
    
    hits <- findOverlaps(Betas, ldmrs, ignore.strand = TRUE)
    ldmr <- Betas[queryHits(hits)]
    mcols(ldmrs) <- cbind(mcols(ldmrs), aggregate(x = mcols(ldmr), by = list(subjectHits(hits)), FUN = mean, na.rm = TRUE)[,-1])
    rm(ldmr, hits)
    
    png(filename = paste0(config$output$outfile_prefix, "_", coef, "_DMR_heatmap.png"), width = 8, height = 11, units = "in", res = 300)

    col_cols <- lapply(anno_columns, function(i) {
        levels <- unique(sample_info[, i])
        coef_cols <- sapply(1:length(levels), function(j) {grey.colors(length(levels))[j]})
        names(coef_cols) <- levels
        coef_cols
    })
    names(col_cols) <- anno_columns
    
    col_anno <- HeatmapAnnotation(df = design, col = col_cols)
    
    print(Heatmap(
        matrix = as.matrix(mcols(ldmrs)[sample_info[ ,"sample"]]),
        cluster_rows = TRUE,
        clustering_distance_rows = "euclidean",
        row_title = NULL,
        cluster_columns = FALSE,
        use_raster = TRUE,
        raster_device = "png",
        top_annotation = col_anno
    ))
    
    dev.off()
    
}

pseudoMAplot <- function(all_cpg, dmrs, coverage, diff, coef, pval_threshold = 1e-5) {
    
    pval_col <- grep(pattern = "pval", names(all_cpg))
    sig <- c(all_cpg[, ..pval_col] <= pval_threshold)
    
    ldmls <- GRanges(
        seqnames = Rle(all_cpg$chr),
        range = IRanges(start = all_cpg$pos, end = all_cpg$pos+1),
        col = c("black", "red")[sig +1],
        pch = c(".", "*")[sig +1]
    )
    
    ldmls <- sort(ldmls)
    
    ldmrs <- GRanges(
        seqnames = Rle(dmrs$chr),
        range = IRanges(start = dmrs$start, end = dmrs$end)
    )
    
    ldmrs <- sort(ldmrs)
    
    hits <- findOverlaps(ldmls, ldmrs, ignore.strand = TRUE)
    
    mcols(ldmls[queryHits(hits)])$col <- "blue"
    
    png(filename = paste0(config$output$outfile_prefix, "_", coef, "_MAplot.png"), width = 8, height = 8, units = "in", res = 300)
    par(mar = c(5, 4, 0, 0) + 0.1)
    plot(x = coverage, 
         y = diff, 
         xlab = "Mean CpG Coverage",
         ylab = "Methylation Difference",
         xlim = c(0, max_cov+5),
         col = "black", 
         pch = ".")
    abline(h = 0, col = "blue")
    points(x = coverage[sig],
           y = diff[sig],
           pch = 16,
           col = mcols(ldmls[sig])$col, 
           cex = 0.3)
    dev.off()
    
}

## Differential Methylation summary and heatmaps

Beta_values <- getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")
ME <- bs_obj_all@rowRanges
mcols(ME) <- as.data.frame(Beta_values)
# mcols(ME) <- as.data.frame(M_values)
rm(Beta_values)

mean_cov <- rowMeans(getCoverage(bs_obj_all, type = "Cov"), na.rm = TRUE)
# var_cov <- rowVars(getCoverage(bs_obj_all, type = "Cov"))

bs_obj_path <- paste0(config$output$outfile_prefix, "_min", min_cov, "_max", max_cov)

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
    
    ## MA plot
    g1 <- levels(factor(samples[, formula_parts]))[1]
    g2 <- levels(factor(samples[, formula_parts]))[2]
    g1 <- samples[samples[, formula_parts] == g1, "sample"]
    g2 <- samples[samples[, formula_parts] == g2, "sample"]
    cov_diff <- rowMeans(as.matrix(mcols(ME)[,g2]), na.rm = TRUE) - rowMeans(as.matrix(mcols(ME)[,g1]), na.rm = TRUE)
    all_cpg <- fread(paste0(bs_obj_path, "_all_sites.txt.gz"))
    dmrs <- fread(paste0(bs_obj_path, "_dmr_delta", delta, "_pval", pval,".txt.gz"))
    pseudoMAplot(all_cpg = all_cpg, dmrs = dmrs, coverage = mean_cov, diff = cov_diff, coef = formula_parts, pval_threshold = pval)
    rm(cov_diff, all_cpg)
    
    dmls <- fread(paste0(bs_obj_path, "_dml_delta", delta, "_pval", pval,".txt.gz"))
    dml_dmr_summary(dmls, dmrs, coef = formula_parts, flag = 0)
    rm(dmls)
    DMR_heatmap(dmrs = dmrs, Betas = ME, design = design, anno_columns = formula_parts, sample_info = samples, coef = formula_parts)

} else if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
    
    for (coef in attr(terms.formula(formula), "term.labels")) {
        if (grepl(":", coef)) {
            coef2 <- strsplit(coef, ":")[[1]]
            coef2 <- lapply(coef2, function(i) {paste0(i, levels(factor(samples[, i]))[2])})
            coef2 <- paste(unlist(coef2), collapse = ".")
        } else {
            coef2 <- paste0(coef, levels(factor(samples[, coef]))[2])
        }
        
        ## MA plot
        g1 <- levels(factor(samples[, coef]))[1]
        g2 <- levels(factor(samples[, coef]))[2]
        g1 <- samples[samples[,coef] == g1, "sample"]
        g2 <- samples[samples[,coef] == g2, "sample"]
        cov_diff <- rowMeans(as.matrix(mcols(ME)[,g2]), na.rm = TRUE) - rowMeans(as.matrix(mcols(ME)[,g1]), na.rm = TRUE)
        all_cpg <- fread(paste0(bs_obj_path, "_", coef2, "_all_sites.txt.gz"))
        dmrs <- fread(paste0(bs_obj_path, "_", coef2, "_dmr_pval", pval,".txt.gz"))
        pseudoMAplot(all_cpg = all_cpg, dmrs = dmrs, coverage = mean_cov, diff = cov_diff, coef = coef2, pval_threshold = pval)
        rm(cov_diff, all_cpg)
        
        dmls <- fread(paste0(bs_obj_path, "_", coef2, "_dml_pval", pval,".txt.gz"))
        dml_dmr_summary(dmls, dmrs, coef = coef2, flag = 1)
        rm(dmls)
        DMR_heatmap(dmrs = dmrs, Betas = ME, design = design, anno_columns = formula_parts, sample_info = samples, coef = coef)
        
    }
}

