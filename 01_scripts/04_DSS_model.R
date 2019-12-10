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
# args <- "~/Projects/safo_epi/methylUtil/config_unpaired.yml"
# args <- "~/Projects/sasa_epi/methylUtil/config_8x8_glm.yml"

## Sanity checking
if (length(args) != 1)
    stop("Usage: DSS_model.R <config.yml>")

if (!is.yaml.file(args[1]))
    stop("You must supply a configuration file in YAML format.\nUsage: DSS_DML_DMR.R <config.yml>")

config <- read.config(args[1])

if (!grepl(config$options$analysis_type, "wald|glm", ignore.case = TRUE))
    stop("Invalid analysis type.")


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
}

if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
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
} else {
    if (config$options$n_cores == 0)
        warning("Using all cores will require a lot of memory")
    n_cores <- ifelse(config$options$n_cores == 0, detectCores(), config$options$n_cores)
    setDTthreads(1)
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

# Fit models
dml_list <- mclapply(chrs, mc.cores = n_cores, mc.preschedule = FALSE, function(chr) {
    
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
    bs_obj <- bs_obj[rowSums(pass) >= min_ind,]
    rm(pass)
    
    # Run linear models
    if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
        # Standard beta-binomial two group test
        message(paste0("Processing chromosome: ", chr))
        capture.output(dml_test <- DMLtest(bs_obj, group1 = grp1, group2 = grp2, smoothing = TRUE))
    }
    
    if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
        # Linear model with family nested in treatment
        message(paste0("Processing chromosome: ", chr))
        capture.output(dml_test <- DMLfit.multiFactor(bs_obj, design = design, formula = formula, smoothing = TRUE))
    }
    
    return(dml_test)
})


## Compile chromosomes, correct global FDR, call DML/DMR, and write results.

# Wald tests
if (grepl(config$options$analysis_type, "wald", ignore.case = TRUE)) {
    dml_test <- do.call(rbind, dml_list)
    dml_test$fdr <- p.adjust(dml_test$pval, method = "BH")
    
    # Write complete outfile...
    fwrite(dml_test, file = paste0(config$output$outfile_prefix, "_all_sites.txt.gz"), quote = FALSE, sep = "\t")
    
    # Call DML and DMR
    dml <- callDML(dml_test, delta = delta, p.threshold = fdr)
    dmr <- callDMR(dml_test, delta = delta, p.threshold = fdr)
    
    # Write DML/DMR outfiles...
    fwrite(dml, file = paste0(config$output$outfile_prefix, "_dml_delta", delta, "_fdr", fdr,".txt.gz"), quote = FALSE, sep = "\t")
    fwrite(dmr, file = paste0(config$output$outfile_prefix, "_dmr_delta", delta, "_fdr", fdr,".txt.gz"), quote = FALSE, sep = "\t")
}

# glm
if (grepl(config$options$analysis_type, "glm", ignore.case = TRUE)) {
    model_mat <- model.matrix(formula, design)
    for (coef in colnames(model_mat)[-1]) {
        dml_factor_test <- mclapply(dml_list, mc.cores = n_cores, mc.preschedule = FALSE, function(chr) {
            test <- DMLtest.multiFactor(chr, coef)
            return(test)
        })
        coef2 <- gsub(":", ".", coef)
        dml_factor_test <- do.call(rbind, dml_factor_test)
        dml_factor_test$fdrs <- p.adjust(dml_factor_test$pval, method = "BH")
        # Write complete outfile...
        fwrite(dml_factor_test, file = paste0(config$output$outfile_prefix, "_", coef2, "_all_sites.txt.gz"), quote = FALSE, sep = "\t")
        
        # Call DML and DMR
        dml <- callDML(dml_factor_test, delta = 0, p.threshold = fdr)
        dmr <- callDMR(dml_factor_test, delta = 0, p.threshold = fdr)
        
        # Write DML/DMR outfiles...
        fwrite(dml, file = paste0(config$output$outfile_prefix, "_", coef2, "_dml_fdr", fdr,".txt.gz"), quote = FALSE, sep = "\t")
        fwrite(dmr, file = paste0(config$output$outfile_prefix, "_", coef2, "_dmr_fdr", fdr,".txt.gz"), quote = FALSE, sep = "\t")
        
    }
}
