#!/usr/bin/env Rscript
## Methylation plotting

for (p in c("configr","tidyverse", "data.table", "BiocManager", "ggplot2", "ggrepel", "ComplexHeatmap", "DSS", "bsseq", "vegan", "parallel")) {
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
# args <- "~/Projects/sasa_epi/methylUtil/config_7x7.yml" ; setwd("~/Projects/sasa_epi/methylUtil")

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
chrs <- read.table(config$input$chrs, header = FALSE, stringsAsFactors = FALSE)[,1]

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

if (is.null(config$options$fdr) | !is.numeric(config$options$fdr)) {
    warning("Invalid FDR. Default to using a value of 0.05.")
    fdr <- 0.05
} else {
    fdr <- config$options$fdr
}

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE) & (is.null(config$options$delta) | !is.numeric(config$options$delta))) {
    warning("Delta required for Wald tests. Default to using a value of 0.1.")
    delta <- 0.1
} else {
    delta <- config$options$delta
}

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
    
    # Filter CpGs on min and max coverage in min individuals
    pass <- getCoverage(bs_obj, type = "Cov") <= max_cov & getCoverage(bs_obj, type = "Cov") >= min_cov
    bs_obj <- bs_obj[rowSums(pass) >= min_ind,]
    rm(pass)
    
    return(bs_obj)
})

bs_obj_all <- suppressWarnings(do.call(rbind, bs_obj_all))


#### MDS plot of samples ####

un <- getCoverage(bs_obj_all, type = "Cov") - getCoverage(bs_obj_all, type = "M")
me <- getCoverage(bs_obj_all, type = "M")

M_values <- log2(me + 1) - log2(un + 1)
colnames(M_values) <- samples[, "sample"]
rm(un, me)

mds <- cmdscale(dist(t(M_values)), k = ncol(M_values) - 1)
colnames(mds) <- paste0("PC",1:ncol(mds))
mds <- data.frame(Sample = row.names(mds), Group = samples[, "group"], mds)
fwrite(mds, "06_methylation_results/methylation_MDS.txt", quote = FALSE, sep = "\t")

all_mds <- ggplot(mds, aes(x = PC1, y = PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_vline(xintercept = 0, colour = "grey") +
    kyles_theme_adjustments +
    geom_text_repel(aes(label = Sample), nudge_y = -10)
ggsave("06_methylation_results/global_methylation_mds_labs.png", plot = all_mds, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)

all_mds_points <- ggplot(mds, aes(x = PC1, y = PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_vline(xintercept = 0, colour = "grey") +
    geom_point(size = 2) +
    kyles_theme_adjustments
ggsave("06_methylation_results/global_methylation_mds_points.png", plot = all_mds_points, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)

rm(mds, all_mds, all_mds_points, M_values)


## Methylation ratios by sample

Beta_values <- getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")
# Beta_values <- asin(sqrt(getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")))

df <- pivot_longer(as.data.frame(Beta_values), cols = starts_with("LJ"), names_to = "Sample", values_to = "Methylation")
rm(Beta_values)
Beta_summary <- df %>% 
    group_by() %>% 
    summarise(Mean = mean(, na.rm = TRUE), Median = median(, na.rm = TRUE), SD = sd(, na.rm = TRUE))

fwrite(Beta_summary, "06_methylation_results/Beta_summary_by_individual.txt", quote = FALSE, sep = "\t")

methyl_ratio_hist <- ggplot(df, aes(x = Methylation)) +
    kyles_theme_adjustments +
    theme(axis.text = element_text(size = 8)) +
    geom_histogram(fill = "steelblue3", binwidth = 0.1) +
    facet_wrap(~Sample, ncol = 4) +
    ylab("Count")
ggsave("03_results/methylation_ratio_hist_by_sample.png", plot = methyl_ratio_hist, device = "png",
       width = 8, height = 8, units = "in", dpi = 300)
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
    
    fwrite(dml_summary, paste0("06_methylation_results/DML_hyper_hypo_distribution", coef, ".txt"), quote = FALSE, sep = "\t")
    
    dmr_summary <- dmrs[, .("areaStat_Mean" = mean(areaStat), 
                            "areaStat_Median" = median(areaStat),
                            "areaStat_SD" = sd(areaStat),
                            "nCG" = sum(nCG),
                            "nCG_mean" = mean(nCG),
                            "nCG_median" = median(nCG),
                            "nCG_SD" = sd(nCG),
                            "N" = .N), by = sign(areaStat)]
    
    fwrite(dml_summary, paste0("06_methylation_results/DMR_hyper_hypo_distribution", coef, ".txt"), quote = FALSE, sep = "\t")
    
}

DMR_heatmap <- function(dmrs, Betas, sample_info, coef = NULL) {
    
    ldmrs <- GRanges(
        seqnames = Rle(dmrs$chr),
        range = IRanges(start = dmrs$start, end = dmrs$end),
        nCG = dmrs$nCG,
        diffMethyl = dmrs$diff.Methy,
        areaStat = dmrs$areaStat
    )
    
    ldmrs <- sort(ldmrs)
    
    hits <- findOverlaps(Betas, ldmrs, ignore.strand = TRUE)
    ldmr <- Betas[queryHits(hits)]
    mcols(ldmrs) <- cbind(mcols(ldmrs), aggregate(x = mcols(ldmr), by = list(subjectHits(hits)), FUN = mean, na.rm = TRUE)[,-1])
    rm(ldmr, hits)
    
    png(filename = paste0("06_methylation_results/", coef, "DMR_heatmap.png"), width = 8, height = 11, units = "in", res = 300)
    Heatmap(
        matrix = as.matrix(mcols(dmrs)[samples_info[,"sample"]]),
        cluster_rows = TRUE,
        #row_split = factor(paste0(seqnames(dmr$dmrs), ":", start(dmr$dmrs), "-", end(dmr$dmrs))),
        row_title = NULL,
        #cluster_row_slices = TRUE,
        cluster_columns = FALSE,
        column_split = as.character(sample_info[,coef]),
        use_raster = TRUE,
        raster_device = "png"
    )
    dev.off()
    
}

## Differential Methylation
if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
    dmls <- fread(paste0(config$output$outfile_prefix, "_dml_delta", delta, "_fdr", fdr,".txt.gz"))
    dmrs <- fread(paste0(config$output$outfile_prefix, "_dmr_delta", delta, "_fdr", fdr,".txt.gz"))
    dml_dmr_summary(dmls, dmrs, flag = 0)
} else if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
    for (coef in as.character(formula)) {
        dmls <- fread(paste0(config$output$outfile_prefix, "_", coef, "_dml_fdr", fdr,".txt.gz"))
        dmrs <- fread(paste0(config$output$outfile_prefix, "_", coef, "_dmr_fdr", fdr,".txt.gz"))
        dml_dmr_summary(dmls, dmrs, coef = coef, flag = 0)
    }
}

## Heatmap

Beta_values <- getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")
ME <- bs_obj_all@rowRanges
mcols(ME) <- as.data.frame(Beta_values)
# mcols(ME) <- as.data.frame(M_values)
rm(Beta_values)




