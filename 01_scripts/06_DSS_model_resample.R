#!/usr/bin/env Rscript
## Working script for DML / DMR quantification using DSS

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "DSS", "bsseq", "parallel")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("DSS", "bsseq")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(library(p, character.only = T))}
    rm(p)
}


args <- commandArgs(T)
# args <- c("~/Projects/sasa_epi/02_data/chrs.txt", "~/Projects/sasa_epi/sample_info.txt", "~/Projects/sasa_epi/03_results/DSS_two_group_test")

if (length(args) != 3)
    stop("Usage: DSS_DML_DMR.R <chr file> <sample info file> <output prefix>")

samples <- read.table(args[2], header = T, stringsAsFactors = FALSE)

if(!identical(colnames(samples), c("sample", "group", "family", "file")))
    stop("Samples file must contain a header row with names: \'sample\', \'group\', \'family\', \'file\'")

design <- data.frame(group = factor(samples[, c('group')], levels = c("Wild", "SAS")), family = samples[, c('family')])


# data_list <- lapply(1:nrow(samples), function(i) {
#     file <- fread(samples[i, "file"], header = FALSE)[,c(-3:-4)]
#     file <- file[V1 %like% c("NC_.*")]
#     file[, V7 := V5 + V6]
#     return(file[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)])
# })
# 
# bs_obj <- makeBSseqData(data_list, samples[,"sample"])
# 
# rm(data_list)
# 

# Simple two group model

chrs <- read.table(args[1], header = FALSE, stringsAsFactors = FALSE)[,1]

n_cores = 8

n_reps = 100

max_cov <- 20
min_cov <- 5

resample <- lapply(1:n_reps, function(i){sample(samples$sample)})

dml_resample_list <- lapply(chrs, function(chr) {
    
    message(paste("Processing chromsome:", chr))

    lsamples <- samples
    
    lsamples$file <- sub("\\.bedGraph\\.gz", paste0("_", chr, "\\.bedGraph\\.gz"), lsamples$file)
    
    data_list <- lapply(1:nrow(samples), function(i) {
        file <- fread(lsamples[i, "file"], header = FALSE)[,c(-3:-4)]
        file[, V7 := V5 + V6]
        return(file[, .("chr" = V1, "pos" = V2, "N" = V7, "X" = V5)])
    })
    
    bs_obj <- makeBSseqData(data_list, samples[,"sample"])
    
    rm(data_list)
    
    #### insert filter routine here
    
    # Min and Max coverage
    pass <- bs_obj@assays$data$Cov <= max_cov & bs_obj@assays$data$Cov >= min_cov
    
    bs_obj <- bs_obj[rowSums(pass) >= 12,]
    
    rm(pass)
    
    # Standard beta-binomial two group test
    dml_list <- mclapply(resample, mc.cores = n_cores, function(i){
        grp1 = i[1:8]
        grp2 = i[9:16]
        dms <- suppressMessages(DMLtest(bs_obj, group1 = grp1, group2 = grp2, smoothing = TRUE))
        dmls <- callDML(dms, delta = 0.2, p.threshold = 0.01)
        dmrs <- suppressWarnings(callDMR(dms, delta = 0.2, p.threshold = 0.01))
        return(list(nloci = ifelse(is.null(nrow(dmls)), 0, nrow(dmls)), nregion = ifelse(is.null(nrow(dmrs)), 0, nrow(dmrs))))
    })
    
    # Linear model with family nested in treatment
    # dml_test_sm <- DMLfit.multiFactor(bs_obj, design = design, formula = ~ group + group:family, smoothing = TRUE) # group + group:family
    dml_list <- as.data.frame(do.call(rbind, dml_list))
    dml_list <- cbind(chr = chr, rep = 1:100, dml_list)
    
    return(dml_list)
})
 
dml_resample <- do.call(rbind, dml_resample_list)

fwrite(dml_resample, file = paste0(args[3], ".txt.gz"), quote = FALSE, sep = "\t")
