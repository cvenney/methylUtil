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
# bs_obj_path <- paste0(config$output$outfile_prefix, "_min", min_cov, "_max", max_cov)

dmr_obj_path <- paste0(config$output$outfile_prefix, "_", names(design), "_DMR_heatmap_mean_Betas.txt")

#### START WGCNA ####
dmrs <- fread(dmr_obj_path)
gr <- GRanges(seqnames = dmrs$seqnames, ranges = IRanges(start = dmrs$start, end = dmrs$end))

## Get arcsin sqrt transformed beta matrix
beta_trans <- as.matrix(asin(sqrt(dmrs[,-1:-5])))

soft_thresh <- pickSoftThreshold(t(beta_trans), networkType = "signed hybrid", corFnc = "bicor", corOptions = list(maxPOutliers = 0.1), powerVector = seq(1,30, by = 1), verbose = 5)

ggplot(soft_thresh$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
    geom_text(aes(label = Power), show.legend = FALSE)
    
ggplot(soft_thresh$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = Power), show.legend = FALSE)

softPower <- 17

adj <- adjacency(t(beta_trans), power = softPower, type = "signed hybrid", corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.1))
dissTOM <- 1 - TOMsimilarity(adj)

tree <- hclust(as.dist(dissTOM), method = "average")

plot(tree, xlab = "", sub = "", main = "CpG clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

minModuleSize <- 3

dynamicMods <- cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)
length(table(dynamicColors))

MEList <- moduleEigengenes(t(beta_trans), colors = dynamicColors)

MEDiss = 1-cor(MEList$eigengenes)

METree = hclust(as.dist(MEDiss), method = "average")

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.2

abline(h = MEDissThres, col = 'red')

merge <- mergeCloseModules(t(beta_trans), dynamicColors, cutHeight = MEDissThres, verbose = 3)

moduleColors <- merge$colors
length(table(moduleColors))

plotDendroAndColors(tree, cbind(dynamicColors, moduleColors),c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

phenos <- read.table("LJ_phenos.txt", header = TRUE, stringsAsFactors = FALSE)

MEs <- orderMEs(moduleEigengenes(t(beta_trans), moduleColors)$eigengenes)

moduleTraitCor <- bicor(MEs, phenos[,-1], maxPOutliers = 0.1, use = "p")

moduleTraitPvalue <- bicorAndPvalue(MEs, phenos[,-1], maxPOutliers = 0.1, use = "p")$p

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 3), ")", sep = "")

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
               main = paste("Module-trait relationships"), cex.lab.y = 0.5)

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(bicor(t(beta_trans), MEs, maxPOutliers = 0.1, use = "p"))
MMPvalue = as.data.frame(bicorAndPvalue(t(beta_trans), MEs, maxPOutliers = 0.1, use = "p")$p)
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(bicor(t(beta_trans), phenos[,-1], maxPOutliers = 0.1, use = "p"))
GSPvalue = as.data.frame(bicorAndPvalue(t(beta_trans), phenos[,-1], maxPOutliers = 0.1, use = "p")$p)
names(geneTraitSignificance) = paste("GS.", names(phenos)[-1], sep="")
names(GSPvalue) = paste("p.GS.", names(phenos)[-1], sep="")

var = 4

which(moduleTraitPvalue[,var] <= 0.1)

# weight / length
module = "pink"
column = match(module, modNames)
moduleGenes = moduleColors==module
sum(moduleGenes)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, var]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

test <- gr[which(moduleGenes)]
dmr_genes <- fread("06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001_dmr_context.txt")
dmr_genes[V9 == "gene", GeneID := sub(";.*", "", sub(".*Name=", "", V15))]
dmr_genes <- unique(GRanges(seqnames = dmr_genes$V1, ranges = IRanges(start = dmr_genes$V2, dmr_genes$V3), geneID = dmr_genes$GeneID))

hits <- findOverlaps(test, dmr_genes)

dmr_genes[subjectHits(hits)]


modTOM = TOMsimilarity(adj)[moduleGenes, moduleGenes]
geneNames <- dmr_genes$geneID[moduleGenes]
geneNames[is.na(geneNames)] <- paste0("NA", 1:sum(is.na(geneNames)))
dimnames(modTOM) = list(geneNames, geneNames)
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = NULL)
## Code to be used with lab-epigenetics script to plot test diffs with samples I actually used
x <- phenos
x <- x %>% mutate(Crosstype = sub("(LJ17)(W|S)(.*)", "\\2", ID), FamilyID = sub("(LJ17)([W|S][0-9])(.*)", "\\2", ID)) %>% 
    mutate(FamilyID = case_when(
        str_starts(FamilyID, "W") ~ paste0("NOR", sub("(W|S)([0-9])", "\\2", FamilyID)), 
        str_starts(FamilyID, "S") ~ paste0("HOR", sub("(W|S)([0-9])", "\\2", FamilyID))))
x$MEpink <- MEs$MEpink
x

ggplot(x, aes(x = MEpink, y = CF, colour = Crosstype)) +
    geom_point()

ggplot(x, aes(x = MEpink, y = CF, colour = Crosstype)) +
    geom_point()
ggplot(x, aes(x = FL, y = WT, colour = Crosstype)) +
    geom_point()
