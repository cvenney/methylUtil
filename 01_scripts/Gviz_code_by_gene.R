#!/usr/bin/env Rscript

# setwd("~/Projects/sasa_epi/methylUtil/")

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "Gviz", "GenomicFeatures", "DSS")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if (p %in% c("DSS", "bsseq", "GenomicRanges")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

options(ucscChromosomeNames = FALSE)

slop <- 10000

## Read in DMRs
adult_dmrs <- fread("06_methylation_results/adults_6x6_min5_max20_groupWild_dmr_pval0.001.txt.gz")
adult_dmrs <- GRanges(seqnames = adult_dmrs$chr, ranges = IRanges(start = adult_dmrs$start, end = adult_dmrs$end))

juvenile_dmrs <- fread("06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001.txt.gz")
juvenile_dmrs <- GRanges(seqnames = juvenile_dmrs$chr, ranges = IRanges(start = juvenile_dmrs$start, end = juvenile_dmrs$end))

genes <- fread("05_bed_files/genes.bed")
genes <- GRanges(seqnames = genes$V1, ranges = IRanges(start = genes$V2, end = genes$V3), strand = genes$V6, geneid = genes$V4)

a_hits <- findOverlaps(adult_dmrs, genes, maxgap = 5000)
j_hits <- findOverlaps(juvenile_dmrs, genes, maxgap = 5000)

a_hits <- data.table(a_hits = queryHits(a_hits), geneid = genes[subjectHits(a_hits)]$geneid)
j_hits <- data.table(j_hits = queryHits(j_hits), geneid = genes[subjectHits(j_hits)]$geneid)

dmr_hits <- merge(a_hits, j_hits, by = "geneid")

juvenile_dmrs <- juvenile_dmrs[dmr_hits$j_hits]
adult_dmrs <- adult_dmrs[dmr_hits$a_hits]

keepGR <- genes[genes$geneid %in% unique(dmr_hits$geneid)]
start(keepGR) <- ifelse(start(keepGR)-slop < 0, 0, start(keepGR)-slop)
end(keepGR) <- end(keepGR)+slop

## Read in and subset raw data
adult_data <- readRDS("06_methylation_results/adults_6x6_min5_max20")
adult_data <- adult_data[subjectHits(findOverlaps(keepGR, adult_data))]
gc()

juvenile_data <- readRDS("06_methylation_results/juveniles_8x8_min5_max20")
juvenile_data <- juvenile_data[subjectHits(findOverlaps(keepGR, juvenile_data))]
gc()

# Build gene models from NCBI annotation
TxDb.Ssalar <- makeTxDbFromGFF("02_reference/genes_with_utrs.gff") 
# TxDb.Ssalar <- makeTxDbFromEnsembl(organism = "Salmo salar", release = 100)

## reorder keepGR for logical order in plot
seqlevels(keepGR) <- sort(seqlevels(keepGR))
keepGR <- sort(keepGR)


pdf("FigureS4_overlapping-DMR.pdf", width = 8, height = 18)
grid.newpage()
vp_y_coord <- 1 + (1 / (length(keepGR)*2))

for (i in 1:length(keepGR)) {
    
    vp_y_coord <- vp_y_coord - (1/ length(keepGR))
    
    pushViewport(viewport(x = 0.5, y = vp_y_coord, height = 1 / length(keepGR), width = 1)) # height = 1/ length(keepGR)

    chr <- as.character(seqnames(keepGR[i]))
    a_sub <- subjectHits(findOverlaps(keepGR[i], adult_dmrs, maxgap = 5000))
    j_sub <- subjectHits(findOverlaps(keepGR[i], juvenile_dmrs, maxgap = 5000))

    genomeAxis <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE,  labelPos = "above", col = "black", fontcolor = "black", fontsize = 6)
    customFromTxDb <- GeneRegionTrack(TxDb.Ssalar,
                                      chromosome = chr,
                                      transcriptAnnotation = "gene",
                                      collapseTranscripts = FALSE,
                                      name = "Genes",
                                      showId = TRUE,
                                      geneSymbol = TRUE,
                                      background.title = "firebrick4",
                                      rotation.title = 0,
                                      shape = c("smallArrow"))

    # Grab individual level methylation data
    ahits <- findOverlaps(keepGR, adult_data)
    abeta <- data.frame(HOR = rowSums(getCoverage(adult_data[subjectHits(ahits), 1:6], type = "M")) / rowSums(getCoverage(adult_data[subjectHits(ahits), 1:6], type = "Cov")),
                        NOR = rowSums(getCoverage(adult_data[subjectHits(ahits), 7:12], type = "M")) / rowSums(getCoverage(adult_data[subjectHits(ahits), 7:12], type = "Cov")))
    # colnames(abeta) <- colnames(adult_data)
    adf <- granges(adult_data[subjectHits(ahits)])
    mcols(adf) <- abeta
    adataTrack <- DataTrack(adf, name = "Adult",
                            groups = rep(c("SAS", "Wild"), each = ncol(abeta)/2),
                            type = c("a"), col = c("darkgoldenrod1", "blue"), lwd = 2,
                            background.title = "firebrick4", ylim = c(0,1.1), showTitle = FALSE, legend = FALSE)
    admrTrack <- AnnotationTrack(reduce(adult_dmrs[a_sub]), name = "Adult", fill = "grey50", col = "black",
                                 background.title = "firebrick4", rotation.title = 0)
    aht <- HighlightTrack(trackList = list(adataTrack, admrTrack),
                          start = start(adult_dmrs[a_sub]),
                          end = end(adult_dmrs[a_sub]),
                          chromosome = chr, col = "grey30", fill = "grey90")

    jhits <- findOverlaps(keepGR, juvenile_data)
    jbeta <- data.frame(HOR = rowSums(getCoverage(juvenile_data[subjectHits(jhits), 1:8], type = "M")) / rowSums(getCoverage(juvenile_data[subjectHits(jhits), 1:8], type = "Cov")),
                        NOR = rowSums(getCoverage(juvenile_data[subjectHits(jhits), 9:16], type = "M")) / rowSums(getCoverage(juvenile_data[subjectHits(jhits), 9:16], type = "Cov")))
    # colnames(jbeta) <- colnames(juvenile_data)
    jdf <- granges(juvenile_data[subjectHits(jhits)])
    mcols(jdf) <- jbeta
    jdataTrack <- DataTrack(jdf, name = "Juvenile",
                            groups = rep(c("SAS", "Wild"), each = ncol(jbeta)/2),
                            type = c("a"), col = c("darkgoldenrod1", "blue"), lwd = 2,
                            background.title = "firebrick4", lty = 2, ylim = c(0,1.1), showTitle = FALSE, legend = FALSE)
    jdmrTrack <- AnnotationTrack(reduce(juvenile_dmrs[j_sub]), name = "Juv.", fill = "grey50", col = "black", lty = 2,
                                 background.title = "firebrick4", rotation.title = 0)
    jht <- HighlightTrack(trackList = list(jdataTrack, jdmrTrack),
                          start = start(juvenile_dmrs[j_sub]),
                          end = end(juvenile_dmrs[j_sub]),
                          chromosome = chr, col = "grey30", fill = "grey90")

    #plot_name <- paste0("06_methylation_results/DMR_",chr, "_", start(keepGR[i]), "-", end(keepGR[i]), ".pdf")

    #pdf(plot_name, width = 7, height = 4)
    plotTracks(trackList = list(genomeAxis, aht, jht, customFromTxDb),
               chromosome = chr,
               main = chr,
               cex.main = 0.5,
               from = start(keepGR[i]),
               to = end(keepGR[i]),
               add = TRUE
    )
    #dev.off()
    
    pushViewport(viewport(x = 0.0325, y = 0.9, width = 0.05, height = 0.05))
    grid.text(label = LETTERS[i], gp = gpar(fontsize = 16, fontface = "bold"))
    popViewport()
    popViewport()
}

dev.off()

# #### Overlaid
# 
# ahits <- findOverlaps(keepGR, adult_data)
# abeta <- getCoverage(adult_data[subjectHits(ahits)], type = "M") / getCoverage(adult_data[subjectHits(ahits)], type = "Cov")
# colnames(abeta) <- colnames(adult_data)
# adf <- granges(adult_data[subjectHits(ahits)])
# mcols(adf) <- abeta
# adataTrack <- DataTrack(adf, name = "Methylation", 
#                         groups = rep(c("HOR", "NOR"), each = ncol(abeta)/2),
#                         type = c("a"), col = c("steelblue3", "firebrick4"),
#                         background.title = "darkblue")
# 
# jhits <- findOverlaps(keepGR, juvenile_data)
# jbeta <- getCoverage(juvenile_data[subjectHits(jhits)], type = "M") / getCoverage(juvenile_data[subjectHits(jhits)], type = "Cov")
# colnames(jbeta) <- colnames(juvenile_data)
# jdf <- granges(juvenile_data[subjectHits(jhits)])
# mcols(jdf) <- jbeta
# jdataTrack <- DataTrack(jdf, name = "Methylation", 
#                         groups = rep(c("HOR", "NOR"), each = ncol(jbeta)/2),
#                         type = c("a"), col = c("steelblue3", "firebrick4"),
#                         background.title = "darkblue", lty = 2)
# 
# ylims <- extendrange(range(c(values(adataTrack), values(jdataTrack))))
# 
# ot <- OverlayTrack(trackList=list(adataTrack, jdataTrack), name = "Methylation", background.title = "darkblue")
# 
# admrTrack <- AnnotationTrack(adult_dmrs[a_sub], name = "Adult DMR", fill = "grey", col = "black", 
#                              background.title = "darkblue")
# 
# jdmrTrack <- AnnotationTrack(juvenile_dmrs[j_sub], name = "Juvenile DMR", fill = "grey", col = "black", lty = 2,
#                              background.title = "darkblue")
# 
# aht <- HighlightTrack(trackList = list(ot, admrTrack),
#                      start = c(start(adult_dmrs[a_sub])), 
#                      end = c(end(adult_dmrs[a_sub])), 
#                      chromosome = chr)
# 
# jht <- HighlightTrack(trackList = list(ot, jdmrTrack),
#                      start = c(start(juvenile_dmrs[j_sub])), 
#                      end = c(end(juvenile_dmrs[j_sub])), 
#                      chromosome = chr)
# 
# plotTracks(trackList = list(genomeAxis, aht, jht, customFromTxDb),
#            chromosome = chr,
#            from = start(keepGR[i]) + 9000, 
#            to = end(keepGR[i]) - 9000
# ) 
# 
# 
