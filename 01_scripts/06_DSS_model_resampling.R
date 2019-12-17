#!/usr/bin/env Rscript
## Working script for DML / DMR quantification using DSS

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "DSS", "bsseq", "parallel", "configr", "tidyverse")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("DSS", "bsseq")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(require(p, character.only = T))}
    rm(p)
}


args <- commandArgs(T)
# args <- c(10, "~/Projects/safo_epi/methylUtil/config_unpaired.yml"); setwd("~/Projects/safo_epi/methylUtil/")
# args <- c(10, "~/Projects/sasa_epi/methylUtil/config_"); setwd("~/Projects/sasa_epi/methylUtil/")

## Sanity checking
if (length(args) != 2)
    stop("Usage: DSS_model_resampling.R <n_reps> <config.yml>")

if (!is.yaml.file(args[2]))
    stop("You must supply a configuration file in YAML format.\nUsage: DSS_DML_DMR.R <n_reps> <config.yml>")

n_reps <- as.integer(args[1])
config <- read.config(args[2])

if (!grepl(config$options$analysis_type, "wald|glm", ignore.case = TRUE))
    stop("Invalid analysis type.")

# if (is.null(config$options$n_replicates) | !is.numeric(config$options$n_replicates)) {
#     stop("Number of replicates not specified. Do you have the correct config file?")
# } else
#     nreps <- config$options$n_replicates

# Parse formula
if (is.null(config$options$formula) | !grepl("\\~", config$options$formula))
    stop("Invalid formula. You must provide a design formula beginning with a tilde (e.g. \'~ Treatment\')")

formula <- as.formula(config$options$formula)
formula_parts <- unlist(strsplit(config$options$formula, split = c("\\~ |\\~|\\s\\+\\s|\\s\\+|\\+\\s|\\+|\\*|\\s\\*|\\s\\*\\s|\\*\\s|\\:")))[-1]

if (length(formula_parts) > 1 & grepl(config$options$analysis_type, "wald", ignore.case = TRUE))
    stop("You specified a Wald test with more than one factor.\nPlease verify your input.")


# Load files and set up design
samples <- read.table(config$input$sample_info, header = T, stringsAsFactors = FALSE)
chrs <- read.table(config$input$chrs, header = FALSE, stringsAsFactors = FALSE)[,1]

if(!all(c("sample", "file", formula_parts) %in% colnames(samples)))
    stop("Samples file must contain a header row with names: \'sample\', \'file\', and given factor(s).")

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
    grp <- formula_parts
    ref <- config$options$reference_condition
    treat <- config$options$treatment_condition
    if (!(ref %in% with(samples, get(grp))) | !(treat %in% with(samples, get(grp))))
        stop("Specified factor levels not found in factor names.\nPlease verify your input.")
    design <- data.frame(group = factor(samples[, grp], levels = c(ref, treat)))
    grp1 = samples[samples[, grp] == levels(design$group)[1], "sample"]
    grp2 = samples[samples[, grp] == levels(design$group)[2], "sample"]
    grp1s <- lapply(1:n_reps, function(i) {sample(samples[, "sample"], length(grp1), replace = TRUE)})
    grp2s <- lapply(1:n_reps, function(i) {sample(samples[, "sample"], length(grp2), replace = TRUE)})
}

if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
    if (any(!formula_parts %in% colnames(samples)))
        stop("Factors specified in formula design are not present in the ")
    design <- data.frame(samples[, formula_parts])
    design[] <- lapply(design, factor)
    resamples <- lapply(1:n_reps, function(i) {sample(samples[, "sample"], nrow(design), replace = TRUE)})
    model_mat <- 
    coefs <- colnames(model.matrix(formula, design))[-1]
}


# Set number of cores for parallel computing
if (is.null(config$options$n_cores)) {
    warning("\'n_cores\' not specified. Default to using 1 core.")
    n_cores <- 1
}
    
if (config$options$n_cores == 0)
    warning("Using all cores will require a lot of memory")

n_cores <- ifelse(config$options$n_cores == 0, detectCores(), config$options$n_cores)

# Set coverage options for filtering
if (is.null(config$options$max_coverage)) {
    warning("\'max_coverage\' not specified. Default to using a value of 30.")
    max_cov <- 30L
} else
    max_cov <- config$options$max_coverage

if (is.null(config$options$max_coverage)) {
    warning("\'max_coverage\' not specified. Default to using a value of 10.")
    min_cov <- 10L
} else
    min_cov <- config$options$min_coverage

if (is.null(config$options$max_coverage)) {
    warning("\'min_individuals\' not specified. Requiring coverage for all individuals.")
    min_ind <- nrow(samples)
} else
    min_ind <- config$options$min_individuals

if (is.null(config$options$pval_threshold) | !is.numeric(config$options$pval_threshold)) {
    warning("Invalid pval_threshold. Default to using a value of 0.05.")
    pval <- 0.05
} else
    pval <- config$options$pval_threshold

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE) & (is.null(config$options$delta) | !is.numeric(config$options$delta))) {
    warning("Delta required for Wald tests. Default to using a value of 0.1.")
    delta <- 0.1
} else
    delta <- config$options$delta


# Fit models
rdmr_list <- lapply(chrs, function(chr) {
    
    # Create local sample info and adujust filenames for specific chr
    lsamples <- samples
    lsamples$file <- sub("\\.bedGraph\\.gz", paste0("_", chr, "\\.bedGraph\\.gz"), lsamples$file)
    
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
    bs_obj <- bs_obj[rowSums(pass, na.rm = TRUE) >= min_ind,]
    rm(pass)
    
    # Run linear models
    if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
        
        message(paste0("Processing chromosome: ", chr))
        rdmrs <- mclapply(1:n_reps, mc.cores = n_cores, function(i) {
            capture.output(dml_test <- DMLtest(bs_obj[,c(grp1s[[i]], grp2s[[i]])], group1 = 1:length(grp1s[[i]]), group2 = (length(grp1s[[i]])+1):(length(grp2s[[i]]) + length(grp1s[[i]])), smoothing = TRUE))
            rdmrs <- callDMR(dml_test, delta = delta, p.threshold = pval)
            if (is.null(rdmrs))
                rdmrs <- GRanges(seqnames = Rle("NULL", 1), IRanges(start = 0, end = 1))
            else
                rdmrs <- GRanges(seqnames = Rle(rdmrs$chr), range = IRanges(start = rdmrs$start, end = rdmrs$end), nCG = rdmrs$nCG, AS = rdmrs$areaStat)
            return(rdmrs)
        })
        
    }
    
    if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
        
        message(paste0("Processing chromosome: ", chr))
        rdmrs <- mclapply(1:n_reps, mc.cores = n_cores, function(i) {
            capture.output(dml_test <- DMLfit.multiFactor(bs_obj[,resamples[[i]]], design = design, formula = formula, smoothing = TRUE))
            dmr_factor_list <- lapply(coefs, function(coef) {
                test <- DMLtest.multiFactor(dml_test, coef)
                rdmrs <- callDMR(test, delta = 0, p.threshold = pval)            
                if (is.null(rdmrs))
                    rdmrs <- GRanges(seqnames = Rle("NULL", 1), IRanges(start = 0, end = 1))
                else
                    rdmrs <- GRanges(seqnames = Rle(rdmrs$chr), range = IRanges(start = rdmrs$start, end = rdmrs$end), nCG = rdmrs$nCG, AS = rdmrs$areaStat)
                return(rdmrs)
            })
            names(dmr_factor_list) <- coefs
            return(dmr_factor_list)
        })
    }
    
    return(rdmrs)
})

if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
    n_DMR <- rowSums(sapply(rdmr_list, function(i) {
        unlist(lapply(1:length(i), function(j) {
            length(i[[j]][seqnames(i[[j]]) != "NULL"])
        })) 
    }, simplify = "matrix"), na.rm = TRUE)

    fwrite(data.frame(n_DMR = n_DMR), paste0(config$output$outfile_prefix, "_resample_n_DMR.txt.gz"), quote = FALSE, sep = "\t")
    
    dmrs <- fread(paste0(config$output$outfile_prefix, "_dmr_delta", delta, "_pval", pval,".txt.gz")) 
    dmrs <- GRanges(
        seqnames = Rle(dmrs$chr),
        range = IRanges(start = dmrs$start, end = dmrs$end)
    )
    
    DMR_overlap <- sapply(1:n_reps, function(i) {
        ldmr <- lapply(rdmr_list, function(j) {
            j[[i]]
        })
        suppressWarnings(ldmr <- do.call(c, ldmr))
        length(suppressWarnings(findOverlaps(ldmr, dmrs)))
    })
    
    fwrite(data.frame(DMR_overlap = DMR_overlap), paste0(config$output$outfile_prefix, "_resample_DMR_overlap.txt.gz"), quote = FALSE, sep = "\t")
}

if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
        
    n_DMR <- lapply(coefs, function(coef) {
        rowSums(sapply(rdmr_list, function(i) {
            unlist(lapply(1:length(i), function(j) {
                length(i[[j]][[coef]][seqnames(i[[j]][[coef]]) != "NULL"])
            })) 
        }, simplify = "matrix"), na.rm = TRUE)
    })
    names(n_DMR) <- coefs
    fwrite(as.data.frame(n_DMR), paste0(config$output$outfile_prefix, "_resample_n_DMR.txt.gz"), quote = FALSE, sep = "\t")
    
    DMR_overlap <- lapply(coefs, function(coef) {
        if (grepl(":", coef))
            coef2 <- gsub(":", "\\.", coef)
        else
            coef2 <- coef
        
        dmrs <- fread(paste0(config$output$outfile_prefix, "_", coef2, "_dmr_pval", pval,".txt.gz")) 
        dmrs <- GRanges(
            seqnames = Rle(dmrs$chr),
            range = IRanges(start = dmrs$start, end = dmrs$end)
        )
        
        sapply(1:n_reps, function(i) {
            ldmr <- lapply(rdmr_list, function(j) {
                j[[i]][[coef]]
            })
            suppressWarnings(ldmr <- do.call(c, ldmr))
            length(suppressWarnings(findOverlaps(ldmr, dmrs)))
        })
    })
    names(DMR_overlap) <- coefs
    fwrite(as.data.frame(DMR_overlap), paste0(config$output$outfile_prefix, "_resample_DMR_overlap.txt.gz"), quote = FALSE, sep = "\t")
}
