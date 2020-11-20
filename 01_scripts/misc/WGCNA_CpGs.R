#!/usr/bin/env Rscript
## Working script for WGCNA analysis of top variable CpG methylation

args <- commandArgs(T)
# args <- "~/Projects/safo_epi/methylUtil/config_unpaired.yml"; setwd("~/Projects/safo_epi/methylUtil")
# args <- "~/Projects/sasa_epi/methylUtil/config_juvenile_samples_8x8.yml"; setwd("~/Projects/sasa_epi/methylUtil")

## Sanity checking
if (length(args) != 1)
    stop("Usage: DSS_model.R <config.yml>")

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "DSS", "bsseq", "GenomicRanges", "dmrseq", "MethCP", "parallel", "configr", "tidyverse", "WGCNA")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if (p %in% c("DSS", "bsseq", "GenomicRanges")) {
            BiocManager::install(p)
        } else if (p == "MethCP") { 
            devtools::install_github(repo = "kylewellband/MethCP", ref = "bug-fixes", repos = BiocManager::repositories())
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

if (!is.yaml.file(args[1]))
    stop("You must supply a configuration file in YAML format.\nUsage: 03_DSS_model.R <config.yml>")

config <- read.config(args[1])

if (!(config$options$analysis_type %in% c("wald", "glm", "dmrseq", "MethCP-wald", "MethCP-glm")))
    stop("Invalid analysis type. Check your spelling. Case-sensitive options are: wald, glm, dmrseq, MethCP-wald, MethCP-glm")


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


#### START WGCNA ####
gr <- granges(bs_obj)

## Get arcsin sqrt transformed beta matrix
beta_trans <- asin(sqrt(getCoverage(bs_obj, type = "M") / getCoverage(bs_obj, type = "Cov")))
cov <- getCoverage(bs_obj, type = "Cov")
rm(bs_obj)

# decreasing_idx_old <- order(rowVars(beta_trans), decreasing = TRUE)
wt_var <- sapply(1:nrow(beta_trans), function(i) {weightedVar(beta_trans[i,], w = cov[i,])})

# ## Fixed number of genes
# num_genes <- 25000
# idx <- order(wt_var, decreasing = TRUE)
# beta_trans <- beta_trans[idx,][1:num_genes,]
# cov <- cov[idx,][1:num_genes,]

# Variance Threshold
var_threshold <- 0.05
idx <- wt_var > var_threshold

# soft_thresh <- pickSoftThreshold(t(beta_trans[idx,]), weights = t(cov[idx,]), networkType = "signed hybrid", powerVector = seq(1,20), verbose = 5)
soft_thresh <- pickSoftThreshold(t(beta[idx,]), 
                                 networkType = "signed hybrid", 
                                 corFnc = "bicor", 
                                 corOptions = list(use = "p", maxPOutliers = 0.1), 
                                 powerVector = seq(1,20), 
                                 verbose = 5)

ggplot(soft_thresh$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
    geom_text(aes(label = Power), show.legend = FALSE) +
    geom_hline(yintercept = 0.8, color = "red")
    
ggplot(soft_thresh$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = Power), show.legend = FALSE)

softPower <- 18

net <- blockwiseModules(t(beta[idx,]),
                        power = softPower,
                        maxBlockSize = 30000,
                        networkType = "signed hybrid",
                        corType = "bicor",
                        maxPOutliers = 0.1,
                        loadTOM = TRUE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "06_methylation_results/blockwiseTOM")

adj <- adjacency(t(beta[idx,]), 
                 power = softPower, 
                 corFnc = "bicor", 
                 corOptions = list(use = "p", maxPOutliers = 0.1))
# dissTOM <- 1 - TOMsimilarity(adj)
# saveRDS(dissTOM, "dissTOM.rds")
dissTOM <- readRDS("dissTOM.rds")

tree <- hclust(as.dist(dissTOM), method = "average")

plot(tree, xlab = "", sub = "", main = "CpG clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

minModuleSize <- 20

dynamicMods <- cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)
length(table(dynamicColors))

MEList <- moduleEigengenes(t(beta[idx,]), colors = dynamicColors)

MEDiss = 1-cor(MEList$eigengenes)

METree = hclust(as.dist(MEDiss), method = "average")

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.5

abline(h = MEDissThres, col = 'red')

merge <- mergeCloseModules(t(beta[idx,]), dynamicColors, cutHeight = MEDissThres, verbose = 3)

moduleColors <- merge$colors
length(table(moduleColors))

plotDendroAndColors(net$dendrograms[[1]], c(moduleColors), c("Dynamic Tree Cut"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

phenos <- read.table("LJ_phenos.txt", header = TRUE, stringsAsFactors = FALSE)

MEs <- orderMEs(moduleEigengenes(t(beta[idx,]), moduleColors)$eigengenes)
# MEs <- net$MEs
# moduleColors <- net$colors

moduleTraitCor <- bicor(MEs, phenos[,-1], use = "p", maxPOutliers = 0.1)
moduleTraitPvalue <- bicorAndPvalue(MEs, phenos[,-1], use = "p")$p

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(phenos)[-1],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(bicor(t(beta[idx,]), MEs, maxPOutliers = 0.1, use = "p"))
MMPvalue = as.data.frame(bicorAndPvalue(t(beta[idx,]), MEs, maxPOutliers = 0.1, use = "p")$p)
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(bicor(t(beta[idx,]), phenos[,-1], maxPOutliers = 0.1, use = "p"))
GSPvalue = as.data.frame(bicorAndPvalue(t(beta[idx,]), phenos[,-1], maxPOutliers = 0.1, use = "p")$p)
names(geneTraitSignificance) = paste("GS.", names(phenos)[-1], sep="")
names(GSPvalue) = paste("p.GS.", names(phenos)[-1], sep="")



for (TOI in names(phenos[,-1])) {
    df <- as.data.frame(t(sapply(names(which(moduleTraitPvalue[,TOI] < 0.05)), function(i) {
        column <- match(substring(i, 3), modNames)
        moduleGenes <- net$colors == substring(i, 3)
        N <- sum(moduleGenes)
        res <- bicorAndPvalue(abs(geneModuleMembership[moduleGenes, column]), 
                              abs(geneTraitSignificance[moduleGenes, paste0("GS.", TOI)]), 
                              use = "p", maxPOutliers = 0.1)
        return(c("module" = i, "nCpG" = N, "cor" = res$bicor, "p" = res$p))
    })), stringsAsFactors = FALSE)
    df$nCpG <- as.integer(df$nCpG)
    df$cor <- as.numeric(df$cor)
    df$p <- as.numeric(df$p)
    assign(TOI, df)
}

var <- 1

which(moduleTraitPvalue[,var] < 0.05)

# weight / length
module = "purple"
module = "skyblue"
module = "burlywood4.1"
module = "indianred.1"
column = match(module, modNames)
moduleGenes = moduleColors==module
sum(moduleGenes)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, var]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)

gr <- GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, end = regions$end))
test <- gr[which(moduleGenes)]

test <- gr[decreasing_idx %in% which(moduleGenes)]

dmrs <- read.table("06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001.bed", header = FALSE, stringsAsFactors = FALSE)
dmrs <- GRanges(seqnames = dmrs$V1, ranges = IRanges(start = dmrs$V2, end = dmrs$V3))

lapply(WT$module, function(mod) {
    test <- gr[net$colors == substring(mod, 3)]
    hits <- findOverlaps(test, dmrs)
    if (length(hits) > 0) {
        return(hits)
    } else {
        return(NULL)
    }
})

name_to_test <- CF$module[9]
test <- gr[idx][net$colors == substring(name_to_test,3)]

table(test@seqnames)

hits <- findOverlaps(test, dmrs)
hits
cols = rep("black", length(test))
cols[queryHits(hits)] <- "red"

test[queryHits(findOverlaps(test,dmrs))]
dmrs[subjectHits(findOverlaps(test,dmrs))]

test[seqnames(test) == "NC_027315.1"]

resamp <- sapply(1:1000, function(i) length(findOverlaps(gr[sample(1:length(gr), sum(moduleGenes), replace = F)], dmrs)))
hist(resamp)
abline(v = length(hits), col = "red")
