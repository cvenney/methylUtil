#!/usr/bin/env Rscript
## Working script for DML / DMR quantification using DSS

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

args <- commandArgs(T)
# args <- "~/Projects/safo_epi/methylUtil/config_unpaired.yml"
# args <- "~/Projects/sasa_epi/methylUtil/config_7x7.yml" ; setwd("~/Projects/sasa_epi/methylUtil")

## Sanity checking
if (length(args) != 1)
    stop("Usage: DSS_model.R <config.yml>")

if (!is.yaml.file(args[1]))
    stop("You must supply a configuration file in YAML format.\nUsage: DSS_DML_DMR.R <config.yml>")

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

#### RDA bewteen groups ####

un <- getCoverage(bs_obj_all, type = "Cov") - getCoverage(bs_obj_all, type = "M")
me <- getCoverage(bs_obj_all, type = "M")

M_values <- log2(me + 1) - log2(un + 1)
colnames(M_values) <- samples[, "sample"]
rm(un, me)

# Beta_values <- getCoverage(bs_obj_all, type = "M") / getCoverage(bs_obj_all, type = "Cov")

data_rda <- t(as.matrix(M_values)); rm(M_values)
# data_rda <- t(as.matrix(Beta_values))

# Simple imputation of NaN with the mean methylation for each position
for (i in 1:ncol(data_rda)) {
    data_rda[is.nan(data_rda[,i]), i] <- mean(data_rda[,i], na.rm = TRUE)
}

rda <- rda(data_rda ~ ., data = design)
fwrite(RsquareAdj(rda), "06_methylation_results/RDA_Rsquared.txt", sep = "\t")

aov_results <- anova.cca(rda, parallel = 1, permutations = how(nperm=199))
fwrite(aov_results, "03_results/RDA_ANOVA_results.txt", sep = "\t")

rm(data_rda, i)

sc <- scores(rda, choices = c(1:2), scaling = 3)
rm(rda)

sc$sites <- data.frame(sc$sites, design)
fwrite(sc$sites, "06_methylation_results/RDA_individual_coords.txt", quote = FALSE, sep = "\t")
sc$sites$Group <- interaction(samples[,c(formula_parts, formula_parts)])

rda_sample_biplot <- ggplot(as.data.frame(sc$sites), aes(x = RDA1, y = PC1, color = Group)) +
    theme_adjustments +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_point(size = 3) +
    #scale_color_manual(values = "") +
    geom_text_repel(aes(label = rownames(as.data.frame(sc$sites))), nudge_x = 0.1, nudge_y = 0.1, size = 2)
ggsave("06_methylation_results/rda_sample_biplot.png", plot = rda_sample_biplot, device = "png",
       width = 4.5, height = 4, units = "in", dpi = 300)
rm(rda_sample_biplot)

rcov <- MASS::cov.rob(sc$species[, "RDA1"]) # robust covariance
resmaha <- mahalanobis(as.matrix(sc$species[, "RDA1"]), rcov$center, rcov$cov) # test statistic
lambda <- median(resmaha) / qchisq(0.5, df = 1)
sc$species <- as.data.frame(sc$species)
sc$species$pvals <- pchisq(resmaha / lambda, 1, lower.tail = FALSE)
rm(resmaha, lambda, rcov)
sc$species$fdr <- p.adjust(sc$species$pvals, method = "BH")
sc$species$chr <- seqnames(bs_obj_all)
fwrite(sc$species, "06_methylation_results/rda_snp_scores.txt.gz", sep = "\t")
sc$species <- setDT(sc$species)
n_keep <- sum(sc$species$fdr < 0.05) + 10000
sc$species <- sc$species[order(abs(RDA1), decreasing = TRUE)][1:n_keep,]

rda_CpG_biplot <- ggplot(sc$species, aes(x = RDA1, y = PC1, colour = fdr < 0.05)) +
    theme_adjustments +
    geom_point() +
    scale_color_manual(values = c("grey", "lightgreen"))
ggsave("06_methylation_results/rda_CpG_biplot.png", plot = rda_CpG_biplot, device = "png",
       width = 5, height = 4, units = "in", dpi = 300)
rm(rda_CpG_biplot)

rda_manhattan <- ggplot(sc$species, aes(x = 1:nrow(sc$species), y = abs(RDA1), colour = fdr < 0.05)) +
   theme_adjustments +
   geom_point(shape = "+")
ggsave("06_methylation_results/rda_manhattan.png", plot = rda_manhattan, device = "png",
      width = 8, height = 4, units = "in", dpi = 300)
rm(sc, rda_manhattan)

