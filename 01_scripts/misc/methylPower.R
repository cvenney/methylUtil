library(data.table)
library(DSS)
library(dmrseq)
library(ggplot2)
library(tidyverse)

#### Functions
makeDMRGrange <- function(DMRresult) {
    gr <- GRanges(seqnames = Rle(values = DMRresult$chr),
                  ranges = IRanges(start = DMRresult$start, end = DMRresult$end),
                  strand = Rle(values = "*", nrow(DMRresult)),
                  nCG = DMRresult$nCG,
                  areaStat = DMRresult$areaStat)
    sort(gr)
}

plotOneDMR <- function(DMR, BSseq, region, flank = 1000) {
    if(!class(DMR)=="GRanges")
        DMR <- makeDMRGrange(DMR)
    
    chr <- as.character(runValue(seqnames(DMR[region,])))
    sn <- colData(BSseq)$sample
    
    lDMR <- DMR[region,]
    start(lDMR) <- start(lDMR) - flank
    end(lDMR) <- end(lDMR) + flank
    
    hits <- suppressWarnings(findOverlaps(lDMR, BSseq@rowRanges))
    pos <- start(BSseq@rowRanges[subjectHits(hits)])

    cov <- data.frame(pos = pos, getCoverage(BSseq[subjectHits(hits)], type = "Cov"))
    colnames(cov)[-1] <- sn
    cov <- pivot_longer(cov, cols = -pos, names_to = "sample", values_to = "cov")
    
    meth <- data.frame(pos = pos, getCoverage(BSseq[subjectHits(hits)], type = "M"))
    colnames(meth)[-1] <- sn
    meth <- pivot_longer(meth, cols = -pos, names_to = "sample", values_to = "M")
    
    df <- merge(cov, meth, by = c("pos", "sample"))
    df$meth <- NA
    df$meth[df$cov > 0] <- df$M[df$cov > 0] / df$cov[df$cov > 0]
    
    
    df <- merge(df, as.data.frame(colData(BSseq)[,c("sample", "group")]), by = "sample")
    df$ticks <- -0.1
    
    dmr_poly <- data.frame(x = c(start(DMR[region,]), end(DMR[region,])),
                           ymin = c(0, 0),
                           ymax = c(1, 1))
    
    gp <- ggplot() +
        theme_linedraw() +
        theme(panel.grid = element_blank()) +
        geom_ribbon(data = dmr_poly, aes(x = x, ymin = ymin, ymax = ymax), fill = "grey90") +
        geom_point(data = df, aes(x = pos, y = ticks), shape = 124) +
        xlim(c(min(df$pos), max(df$pos))) +
        ylim(c(-0.1,1)) +
        labs(x = chr, y = "Methylation", size = "Coverage", colour = "Temperature") +
        geom_hline(yintercept = c(0,1), colour = "grey30") +
        geom_point(data = df, aes(x = pos, y = meth, size = cov, colour = group), alpha = 0.4) +
        geom_smooth(data = df, aes(x = pos, y = meth, group = group, colour = group), se = TRUE, method = lm, formula = y ~ splines::bs(x, degree = 4)) #+
        #scale_color_brewer(type = "div", palette = 1)
    
    print(gp)
}

#### Aquire analysis options

# args <- commandArgs(TRUE)
args <- c("~/Projects/safo_epi/methylUtil/06_methylation_results/unpaired_control_min5_max30_all_sites.txt.gz")

#### Read in fitted model results for CpGs

dml_test <- fread(args[1])

dmr <- callDMR(dml_test, delta = 0.2, p.threshold = 0.001)

dmr_gr <- makeDMRGrange(dmr)

hits <- findOverlaps(dmr_gr, bs_obj@rowRanges)


#### Thresholds...

par(mfrow = c(1,3))

hist(dmr$nCG)
hist(dmr$diff.Methy)
hist(dmr$areaStat)
