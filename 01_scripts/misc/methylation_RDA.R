#!/usr/bin/env Rscript
## Working script for RDA analysis of WGBS data

args <- commandArgs(T)
# args <- "config_juvenile_samples_8x8.yml" ; setwd("~/Projects/sasa_epi/methylUtil")

## Sanity checking
if (length(args) != 1)
    stop("Usage: DSS_model.R <config.yml>")


## Install and load necessary packages
for (p in c("data.table", "BiocManager", "DSS", "bsseq", "parallel", "configr", "tidyverse", "vegan", "ggplot2", "ggrepel")) {
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

if (!is.yaml.file(args[1]))
    stop("You must supply a configuration file in YAML format.\nUsage: DSS_DML_DMR.R <config.yml>")

config <- read.config(args[1])


# Parse formula
if (is.null(config$options$formula) | !grepl("\\~", config$options$formula))
    stop("Invalid formula. You must provide a design formula beginning with a tilde (e.g. \'~ Treatment\')")

formula <- as.formula(config$options$formula)
formula_parts <- unlist(strsplit(config$options$formula, split = c("\\~ |\\~|\\s\\+\\s|\\s\\+|\\+\\s|\\+|\\*|\\s\\*|\\s\\*\\s|\\*\\s|\\:")))[-1]

if (length(formula_parts) > 1 & grepl(config$options$analysis_type, "wald|MethCP-wald", ignore.case = TRUE))
    stop("You specified a Wald test with more than one factor.\nPlease verify your input.")


# Load files and set up design
samples <- read.table(config$input$sample_info, header = T, stringsAsFactors = FALSE)

if(!all(c("sample", "file", formula_parts) %in% colnames(samples)))
    stop("Samples file must contain a header row with names: \'sample\', \'file\', and given factor(s).")

if (config$options$analysis_type %in% c("wald", "MethCP-wald")) {
    grp <- formula_parts
    ref <- config$options$reference_condition
    treat <- config$options$treatment_condition
    if (!(ref %in% with(samples, get(grp))) | !(treat %in% with(samples, get(grp))))
        stop("Specified factor levels not found in factor names.\nPlease verify your input.")
    design <- data.frame(group = factor(samples[, grp], levels = c(ref, treat)))
    grp1 = samples[samples[, grp] == levels(design$group)[1], "sample"]
    grp2 = samples[samples[, grp] == levels(design$group)[2], "sample"]
}

if (config$options$analysis_type %in% c("glm", "MethCP-glm")) {
    if (any(!formula_parts %in% colnames(samples)))
        stop("Factors specified in formula design are not present in the ")
    design <- data.frame(samples[, formula_parts])
    design[] <- lapply(design, factor)
    names(design) <- formula_parts
}


# Set number of cores for parallel computing
if (is.null(config$options$n_cores)) {
    warning("\'n_cores\' not specified. Default to using 1 core.")
    n_cores <- 1
    setDTthreads(n_cores)
    bp_params <- MulticoreParam(workers = n_cores)
} else {
    if (config$options$n_cores == 0)
        warning("Using all cores will require a lot of memory")
    n_cores <- ifelse(config$options$n_cores == 0, detectCores(), config$options$n_cores)
    setDTthreads(n_cores)
    bp_params <- MulticoreParam(workers = n_cores)
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
    pval <- config$options$pval_threshold
}

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE) & (is.null(config$options$delta) | !is.numeric(config$options$delta))) {
    warning("Delta required for Wald tests. Default to using a value of 0.1.")
    delta <- 0.1
} else {
    delta <- config$options$delta
}


## Load and save BSseq object
bs_obj_path <- paste0(config$output$outfile_prefix, "_min", min_cov, "_max", max_cov)
if (file.exists(bs_obj_path)) {
    message("Loading existing BSseq object")
    bs_obj <- readRDS(file = bs_obj_path)
} else if (file.exists(paste0("06_methylation_results/", gsub("\\..*", "", config$input$sample_info), "_all_data.rds"))) {
    message("Loading existing BSseq object")
    bs_obj_all <- readRDS(file = paste0("06_methylation_results/", gsub("\\..*", "", config$input$sample_info), "_all_data.rds"))
    keep <- (rowSums(getCoverage(bs_obj_all, type = "Cov") >= min_cov & getCoverage(bs_obj_all, type = "Cov") <= max_cov)) >= min_ind
    bs_obj <- bs_obj_all[keep, ]
    rm(bs_obj_all, keep)
    message("Saving BSseq obj for future use...")
    saveRDS(object = bs_obj, file = bs_obj_path, compress = "gzip")
} else {    
    bs_obj_all <- lapply(1:nrow(samples), function(i) {
        message(paste0("Loading sample: ", samples[i, "sample"]))
        samp <- fread(samples[i, "file"], header = FALSE)[,c(-3:-4)]
        samp[, V7 := V5 + V6] # Combine me+ and me- counts for total coverage
        #bs_obj <- makeBSseqData(list(samp[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)]), samples[i,"sample"])
        return(samp[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)])
    })
    bs_obj_all <- suppressWarnings(makeBSseqData(dat = bs_obj_all, sampleNames = samples$sample))
    bs_obj_all <- bs_obj_all[(rowSums(getCoverage(bs_obj_all, type = "Cov") >= 1) == ncol(bs_obj_all)), ]
    saveRDS(bs_obj_all, paste0("06_methylation_results/", gsub("\\..*", "", basename(config$input$sample_info)), "_all_data.rds"), compress = "gzip")
    
    # Filter CpGs on min and max coverage in min individuals
    keep <- (rowSums(getCoverage(bs_obj_all, type = "Cov") >= min_cov & getCoverage(bs_obj_all, type = "Cov") <= max_cov)) >= min_ind
    bs_obj <- bs_obj_all[keep, ]
    rm(bs_obj_all, keep)
    message("Saving BSseq obj for future use...")
    saveRDS(object = bs_obj, file = bs_obj_path, compress = "gzip")
}

if (file.exists("02_reference/chrs.txt")) {
    chr <- fread("02_reference/chrs.txt", header = FALSE)
    chr <- chr[V1 %in% as.character(runValue(seqnames(bs_obj)))] 
} else {
    chr <- data.table(V1 = as.character(runValue(seqnames(bs_obj))))
}


#### RDA bewteen groups ####

bs_obj <- bs_obj[seqnames(bs_obj) %in% chr$V1]
gr <- granges(bs_obj)

un <- getCoverage(bs_obj, type = "Cov") - getCoverage(bs_obj, type = "M")
me <- getCoverage(bs_obj, type = "M")

M_values <- log2(me + 1) - log2(un + 1)
# M_values <- getCoverage(bs_obj, type = "M") / getCoverage(bs_obj, type = "Cov")

colnames(M_values) <- samples[, "sample"]
rm(un, me, bs_obj)

M_values <- as(M_values, "HDF5Matrix")

## Calculate variance by CpG
setAutoGridMaker(GRIDMAKER = "rowGrid")
var <- unlist(blockApply(M_values, rowVars, BPPARAM = MulticoreParam(workers = n_cores)))
keep <- which(var >= quantile(var, 0.8))

gr <- gr[keep]

data_rda <- t(as.matrix(M_values[keep,]))

# # Simple imputation of NaN with the mean methylation for each position
# for (i in 1:ncol(data_rda)) {
#     data_rda[is.nan(data_rda[,i]), i] <- mean(data_rda[,i], na.rm = TRUE)
# }

## distance-based RDA
pcoa <- cmdscale(dist(data_rda, method = "euclidian"), k = nrow(data_rda)-1, eig = TRUE)
eig_var <-pcoa$eig[-c(nrow(data_rda))]/sum(pcoa$eig[-c(nrow(data_rda))])
pco_keep <- eig_var > bstick(n = length(pcoa$eig)-1)
if (!pco_keep[1])
    pco_keep <- cumsum(eig_var) < 0.9

dbrda <- rda(pcoa$points[, pco_keep] ~ ., data = design)
dbrda_aov_results <- anova.cca(dbrda, parallel = 1, permutations = how(nperm=999))
fwrite(RsquareAdj(dbrda), "06_methylation_results/dbRDA_Rsquared.txt", sep = "\t")
fwrite(dbrda_aov_results, "03_results/dbRDA_ANOVA_results.txt", sep = "\t")

dbsc <- scores(dbrda, choices = c(1:2), scaling = 3)
dbsc$sites <- data.frame(dbsc$sites, design)
dbsc$sites$Group <- interaction(samples[,c(formula_parts, formula_parts)])

dbrda_sample_biplot <- ggplot(as.data.frame(dbsc$sites), aes(x = !!sym(colnames(dbsc$sites)[1]), y = !!sym(colnames(dbsc$sites)[2]), colour = Group)) +
    theme_adjustments +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_point(size = 3) +
    #scale_color_manual(values = "") +
    labs(x = paste0(colnames(dbsc$sites)[1], " (", round((eigenvals(dbrda)/sum(eigenvals(dbrda))*100)[1], digits = 2), "%)"),
         y = paste0(colnames(dbsc$sites)[2], " (", round((eigenvals(dbrda)/sum(eigenvals(dbrda))*100)[2], digits = 2), "%)")) +
    geom_text_repel(aes(label = rownames(as.data.frame(dbsc$sites))), nudge_x = 0.1, nudge_y = 0.1, size = 2)
ggsave("06_methylation_results/dbRDA_sample_biplot.png", plot = dbrda_sample_biplot, device = "png",
       width = 4.5, height = 4, units = "in", dpi = 300)
rm(dbrda_sample_biplot)


## RDA
rda <- rda(data_rda ~ ., data = design)
rm(data_rda, i)

fwrite(RsquareAdj(rda), "06_methylation_results/RDA_Rsquared.txt", sep = "\t")

aov_results <- anova.cca(rda, parallel = 1, permutations = how(nperm=199))
fwrite(aov_results, "03_results/RDA_ANOVA_results.txt", sep = "\t")

sc <- scores(rda, choices = c(1:2), scaling = 3)
rm(rda)

sc$sites <- data.frame(sc$sites, design)
fwrite(sc$sites, "06_methylation_results/RDA_individual_coords.txt", quote = FALSE, sep = "\t")
sc$sites$Group <- interaction(samples[,c(formula_parts, formula_parts)])

rda_sample_biplot <- ggplot(as.data.frame(sc$sites), aes(x = !!sym(colnames(sc$sites)[1]), y = !!sym(colnames(sc$sites)[2]), colour = Group)) +
    theme_adjustments +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_point(size = 3) +
    #scale_color_manual(values = "") +
    labs(x = paste0(colnames(sc$sites)[1], " (", round((eigenvals(rda)/sum(eigenvals(rda))*100)[1], digits = 2), "%)"),
         y = paste0(colnames(sc$sites)[2], " (", round((eigenvals(rda)/sum(eigenvals(rda))*100)[2], digits = 2), "%)")) +
    geom_text_repel(aes(label = rownames(as.data.frame(sc$sites))), nudge_x = 0.1, nudge_y = 0.1, size = 2)
ggsave("06_methylation_results/RDA_sample_biplot.png", plot = rda_sample_biplot, device = "png",
       width = 4.5, height = 4, units = "in", dpi = 300)
rm(rda_sample_biplot)

rcov <- MASS::cov.rob(sc$species[, "RDA1"]) # robust covariance
resmaha <- mahalanobis(as.matrix(sc$species[, "RDA1"]), rcov$center, rcov$cov) # test statistic
lambda <- median(resmaha) / qchisq(0.5, df = 1)
sc$species <- as.data.frame(sc$species)
sc$species$pvals <- pchisq(resmaha / lambda, 1, lower.tail = FALSE)
rm(resmaha, lambda, rcov)
sc$species$fdr <- p.adjust(sc$species$pvals, method = "BH")
sc$species$chr <- as.character(seqnames(gr))
sc$species$pos <- as.integer(start(gr))
fwrite(sc$species, "06_methylation_results/RDA_snp_scores.txt.gz", sep = "\t")
sc$species <- setDT(sc$species)
n_keep <- sum(sc$species$fdr < 0.05) + 10000

rda_CpG_biplot <- ggplot(sc$species[order(abs(RDA1), decreasing = TRUE)][1:n_keep,], 
                         aes(x = !!sym(colnames(sc$sites)[1]), y = !!sym(colnames(sc$sites)[2]), colour = fdr < 0.05)) +
    theme_adjustments +
    geom_point() +
    scale_color_manual(values = c("grey", "lightgreen"))
ggsave("06_methylation_results/RDA_CpG_biplot.png", plot = rda_CpG_biplot, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)
rm(rda_CpG_biplot)

rda_manhattan <- ggplot(sc$species, aes(x = 1:nrow(sc$species), y = abs(RDA1), colour = fdr < 0.05)) +
   theme_adjustments +
   geom_point(shape = "+")
ggsave("06_methylation_results/RDA_manhattan.png", plot = rda_manhattan, device = "png",
      width = 8, height = 4, units = "in", dpi = 300)
rm(sc, rda_manhattan)

