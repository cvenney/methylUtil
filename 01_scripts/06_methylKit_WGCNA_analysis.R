#!/usr/bin/env Rscript

for (p in c("data.table", "BiocManager", "methylKit", "WGCNA", "ggplot2", "cowplot", "network", "sna", "GGally")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if (p %in% c("methylKit")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://cloud.r-project.org", dependencies = T)}
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

# setwd("~/Projects/sasa_epi/methylUtil")

files <- list.files("04_filtered_bedGraphs/", pattern = "methylKit", full.names = TRUE)

file_list <- lapply(files, function(i) i)

sample_names <- lapply(files, function(i) sub("(.*/)(.*)(\\..*).methylKit.gz", "\\2", i))

myobjDB <- methRead(file_list,
                 sample.id = sample_names,
                 assembly = "ICSASG_V2",
                 treatment = rep(c(1,0), each = 8),
                 context = "CpG",
                 dbtype = "tabix",
                 dbdir = "10_methylDB",
                 mincov = 1
)

filtered_myobj <- filterByCoverage(myobjDB, lo.count=1, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

meth <- unite(filtered_myobj, destrand = TRUE)

tiles <- tileMethylCounts(meth, win.size = 100, step.size = 100, cov.bases = 3)

# methyl_diff <- calculateDiffMeth(tiles, adjust = "SLIM")

# WGCNA
regions <- getData(tiles)[wt_var > 0.05, 1:3]
gr <- GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, end = regions$end))

cov <- as.matrix(getData(tiles)[,tiles@coverage.index])
beta <- as.matrix(getData(tiles)[,tiles@numCs.index] / cov)

wt_var <- rowVars(beta)
# Variance Threshold
var_threshold <- 0.05
hist(wt_var, breaks = 100, xlim = c(0, (2 * var_threshold)))
abline(v = var_threshold, col = "red")
sum(wt_var > 0.05)

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

MEs <- net$MEs
row.names(MEs) <- tiles@sample.ids
fwrite(MEs, "06_methylation_results/blockwiseMEs.txt.gz", sep = "\t", row.names = TRUE)

moduleColors <- net$colors
fwrite(cbind(regions, moduleColors), "06_methylation_results/blockwiseModuleColors.txt", sep = "\t", row.names = TRUE)

# png("FigureS4_moduledendros.png", width = 8, height = 12, units = "in", res = 300)
# grid.newpage()
# vp_y_coord <- 1 + (1 / (length(net$dendrograms)*2))
# for (i in 1:length(net$dendrograms)) {
#     vp_y_coord <- vp_y_coord - 1/ length(net$dendrograms)
#     pushViewport(viewport(x = 0.5, y = vp_y_coord, height = 1/length(net$dendrograms), width = 1))
#     grid.rect()
#     par("plt" = gridPLT())
#     plotDendroAndColors(net$dendrograms[[i]], colors = moduleColors[net$blockGenes[[i]]], dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, new = FALSE)
#     popViewport()
# }
# dev.off()

phenos <- read.table("LJ_phenos.txt", header = TRUE, stringsAsFactors = FALSE)

moduleTraitCor <- bicor(MEs, phenos[,-1], use = "p", maxPOutliers = 0.1)
moduleTraitPvalue <- bicorAndPvalue(MEs, phenos[,-1], use = "p")$p

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

sig_cors <- which(rowSums(moduleTraitPvalue < 0.05) > 0)

pdf("methylModule_trait_correlation_heatmap.pdf", width = 8, height = 8)
labeledHeatmap(Matrix = moduleTraitCor[sig_cors,],
               xLabels = names(phenos)[-1],
               yLabels = substring(names(MEs), 3)[sig_cors],
               ySymbols = names(MEs)[sig_cors],
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[sort(c(sig_cors, (sig_cors + nrow(moduleTraitCor)), (sig_cors + 2*nrow(moduleTraitCor)), (sig_cors + 3*nrow(moduleTraitCor))))],
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(bicor(t(beta[idx,]), MEs, maxPOutliers = 0.1, use = "p"))
MMPvalue = as.data.frame(bicorAndPvalue(t(beta[idx,]), MEs, maxPOutliers = 0.1, use = "p")$p)
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
fwrite(cbind(regions, geneModuleMembership, MMPvalue), "06_methylation_results/blockwiseMM.txt.gz", sep = "\t", row.names = TRUE)

geneTraitSignificance = as.data.frame(bicor(t(beta[idx,]), phenos[,-1], maxPOutliers = 0.1, use = "p"))
GSPvalue = as.data.frame(bicorAndPvalue(t(beta[idx,]), phenos[,-1], maxPOutliers = 0.1, use = "p")$p)
names(geneTraitSignificance) = paste("GS.", names(phenos)[-1], sep="")
names(GSPvalue) = paste("p.GS.", names(phenos)[-1], sep="")
fwrite(cbind(regions, geneTraitSignificance, GSPvalue), "06_methylation_results/blockwiseGS.txt.gz", sep = "\t", row.names = TRUE)

correlated_mods <- lapply(names(phenos[,-1]), function(TOI) {
    df <- as.data.frame(t(sapply(names(which(moduleTraitPvalue[,TOI] < 0.05)), function(i) {
        column <- match(substring(i, 3), modNames)
        moduleGenes <- moduleColors == substring(i, 3)
        N <- sum(moduleGenes)
        res <- bicorAndPvalue(abs(geneModuleMembership[moduleGenes, column]), 
                              abs(geneTraitSignificance[moduleGenes, paste0("GS.", TOI)]), 
                              use = "p", maxPOutliers = 0.1)
        return(c("module" = i, "nCpG" = N, "cor" = res$bicor, "p" = res$p))
    })), stringsAsFactors = FALSE)
    df$nCpG <- as.integer(df$nCpG)
    df$cor <- as.numeric(df$cor)
    df$p <- as.numeric(df$p)
    df
})
names(correlated_mods) <- names(phenos[,-1])

dmrs <- read.table("06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001.bed", header = FALSE, stringsAsFactors = FALSE)
dmrs <- GRanges(seqnames = dmrs$V1, ranges = IRanges(start = dmrs$V2, end = dmrs$V3))
genes <- fread("05_bed_files/genes.bed", header = F)
genes <- GRanges(seqnames = genes$V1, ranges = IRanges(start = genes$V2, end = genes$V3), strand = genes$V6, geneid = genes$V4)

correlated_mod_dmrs <- lapply(correlated_mods, function(toi) {
    sapply(toi$module, function(mod) {
        test <- gr[net$colors == substring(mod, 3)]
        hits <- findOverlaps(test, dmrs)
        if (length(hits) > 0) {
            return(hits)
        } else {
            return(NULL)
        }})
})

var <- "CF"
module <- "darkorange"
test <- gr[net$colors == module]

hits <- findOverlaps(test, dmrs)
cols = rep("black", length(test))
cols[queryHits(hits)] <- "red"

ghits <- findOverlaps(test, genes)
chars <- rep(1, length(test))
chars[queryHits(ghits)] <- 19

column = match(module, modNames)
moduleGenes = moduleColors==module
sum(moduleGenes)

pdf(paste0("MM_", module, "_GS_", var, "_regions_DMR_overlap.pdf"), width = 8, height = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, paste0("GS.", var)]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", var),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = cols, pch = chars)
dev.off()


# Resample code to assess significance of DMRs overlapped
sapply(1:1000, function(i){length(unique(subjectHits(findOverlaps(gr[sample(1:length(gr), 100, replace = F)], dmrs))))})



### Plotting etc.
theme_adjustments <- theme_linedraw() + theme(axis.text = element_text(size = 8, colour = "black"), 
                                              axis.title = element_text(size = 10, colour = "black"),
                                              panel.grid = element_blank(), panel.border = element_rect(colour = "black"))
makeGR <- function(df) {
    GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
}

MEs <- fread("06_methylation_results/blockwiseMEs.txt.gz", sep = "\t")
MEs <- as.matrix(MEs[,-1], rownames = MEs$V1)
moduleColors <- fread("06_methylation_results/blockwiseModuleColors.txt", sep = "\t")
regions <- moduleColors[, .(chr, start, end)]
moduleColors <- moduleColors$moduleColors
geneModuleMembership <- as.matrix(fread("06_methylation_results/blockwiseMM.txt.gz", sep = "\t")[, 5:(ncol(MEs)+4)])
geneTraitSignificance <- as.matrix(fread("06_methylation_results/blockwiseGS.txt.gz", sep = "\t")[,5:8])

## Purple WT
purple_idx <- which(moduleColors == "purple")
purple <- data.frame(regions[purple_idx,], 
                     MMpurple = geneModuleMembership[purple_idx, c("MMpurple")], 
                     geneTraitSignificance[purple_idx, c("GS.WT", "GS.FL")], stringsAsFactors = FALSE)
purple$DMR <- "No"
purple$DMR[queryHits(findOverlaps(makeGR(purple), dmrs))] <- "Yes"
plot_purple <- ggplot(purple, aes(x = MMpurple, y = abs(GS.FL))) +
    geom_point(aes(colour = DMR), show.legend = F) + geom_smooth(method = "lm", colour = "black") +
    scale_color_manual(values = c(No = "grey60", Yes = "red")) +
    theme_adjustments + labs(x = "Purple MM", y = "Fork Length GS")
t.test(MMpurple ~ DMR, purple)
t.test(GS.WT ~ DMR, purple)

## Yellow4 CF
y4_idx <- which(moduleColors == "yellow4")
y4 <- data.frame(regions[y4_idx,], 
                     MMyellow4 = geneModuleMembership[y4_idx, c("MMyellow4")], 
                     GS.CF = geneTraitSignificance[y4_idx, c("GS.CF")], stringsAsFactors = FALSE)
y4$DMR <- "No"
y4$DMR[queryHits(findOverlaps(makeGR(y4), dmrs))] <- "Yes"
plot_y4 <- ggplot(y4, aes(x = MMyellow4, y = abs(GS.CF))) +
    geom_point(aes(colour = DMR), show.legend = F) + geom_smooth(method = "lm", colour = "black") +
    scale_color_manual(values = c(No = "grey60", Yes = "red")) +
    theme_adjustments + labs(x = "Yellow4 MM", y = "Condition Factor GS")
t.test(MMyellow4 ~ DMR, y4)
t.test(GS.CF ~ DMR, y4)


## darkorange CF
dor_idx <- which(moduleColors == "darkorange")
dor <- data.frame(regions[dor_idx,], 
                 MMdarkorange = geneModuleMembership[dor_idx, c("MMdarkorange")], 
                 GS.CF = geneTraitSignificance[dor_idx, c("GS.CF")], stringsAsFactors = FALSE)
dor$DMR <- "No"
dor$DMR[queryHits(findOverlaps(makeGR(dor), dmrs))] <- "Yes"
plot_dor <- ggplot(dor, aes(x = MMdarkorange, y = abs(GS.CF))) +
    geom_point(aes(colour = DMR)) + geom_smooth(method = "lm", colour = "black") +
    scale_color_manual(values = c(No = "grey60", Yes = "red")) +
    theme_adjustments + labs(x = "Dark orange MM", y = "Condition Factor GS")
t.test(MMdarkorange ~ DMR, dor)
t.test(GS.CF ~ DMR, dor)


## navajowhite1 CF
nw1_idx <- which(moduleColors == "navajowhite1")
nw1 <- data.frame(regions[nw1_idx,], 
                  MMnavajowhite1 = geneModuleMembership[nw1_idx, c("MMnavajowhite1")], 
                  GS.CF = geneTraitSignificance[nw1_idx, c("GS.CF")], stringsAsFactors = FALSE)
nw1$Genes <- "No"
nw1$Genes[queryHits(findOverlaps(makeGR(nw1), genes, maxgap = 5000))] <- "Yes"
nw1$GeneNames[queryHits(findOverlaps(makeGR(nw1), genes, maxgap = 5000))] <- genes[subjectHits(findOverlaps(makeGR(nw1), genes, maxgap = 5000))]$geneid
plot_nw1 <- ggplot(nw1, aes(x = MMnavajowhite1, y = abs(GS.CF))) +
    geom_point(aes(colour = Genes)) + geom_smooth(method = "lm", colour = "black") +
    scale_color_manual(values = c(No = "grey60", Yes = "red")) +
    theme_adjustments + labs(x = "Navajowhite1 MM", y = "Condition Factor GS")
t.test(MMnavajowhite1 ~ Genes, nw1)
t.test(GS.CF ~ Genes, nw1)

pg <- plot_grid(plot_purple, plot_y4, plot_dor, nrow = 1, rel_widths = c(1,1,1.3), labels = "AUTO")
save_plot("Figure4_GS_MM_plots.pdf", pg, base_height = 3, base_width = 9)

## Write gene IDs associtaed with modules
for (m in names(sig_cors)) {
    print(m)
    fwrite(list(genes[subjectHits(findOverlaps(makeGR(regions[moduleColors == substring(m,3)]), genes, maxgap = 5000))]$geneid), paste0("geneids_module_", substring(m,3), ".txt"), col.names = F)
}

#### Quick test with intra-modular connectivity

# tiles <- fread("/mnt/HDD2/sasa_epi/10_methylDB/methylBase_723c196439c7_tiled.txt")
names(tiles)[1:3] <- names(regions)
tiles <- merge(regions, tiles)
cov <- as.matrix(tiles[,seq(5, ncol(tiles), by = 3), with = FALSE])
beta <- as.matrix(tiles[,seq(6, ncol(tiles), by = 3), with = FALSE] / cov)
# adj <- adjacency(t(beta), power = 18, type = "signed hybrid", corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.1))


purple_net <- network.initialize(length(purple_idx), directed = FALSE)
purple_adj <- bicor(t(beta[purple_idx,]), use = "p", maxPOutliers = 0.1)
purple_adj[purple_adj < 0.8] <- 0
purple_adj <- purple_adj^12
network.adjacency(purple_adj, purple_net, ignore.eval = FALSE, names.eval = "adj")
purple_net %v% "DMR" <- purple$DMR
plot_purple_net <- ggnet2(purple_net, color = "DMR", edge.size = "adj", size = 3, palette = c("Yes" = "red", "No" = "grey60"), legend.position = "none")

y4_net <- network.initialize(length(y4_idx), directed = FALSE)
y4_adj <- bicor(t(beta[y4_idx,]), use = "p", maxPOutliers = 0.1)
y4_adj[y4_adj < 0.8] <- 0
y4_adj <- y4_adj^12
network.adjacency(y4_adj, y4_net, ignore.eval = FALSE, names.eval = "adj")
y4_net %v% "DMR" <- y4$DMR
plot_y4_net <- ggnet2(y4_net, color = "DMR", edge.size = "adj", size = 3, palette = c("Yes" = "red", "No" = "grey60"), legend.position = "none")

dor_net <- network.initialize(length(dor_idx), directed = FALSE)
dor_adj <- bicor(t(beta[dor_idx,]), use = "p", maxPOutliers = 0.1)
dor_adj[dor_adj < 0.8] <- 0
dor_adj <- dor_adj^12
network.adjacency(dor_adj, dor_net, ignore.eval = FALSE, names.eval = "adj")
dor_net %v% "DMR" <- dor$DMR
plot_dor_net <- ggnet2(dor_net, color = "DMR", edge.size = "adj", size = 3, palette = c("Yes" = "red", "No" = "grey60"), legend.position = "none")


pg <- plot_grid(plot_purple, plot_y4, plot_dor, plot_purple_net, plot_y4_net, plot_dor_net, nrow = 2, rel_widths = c(1,1,1.3), labels = "AUTO")
save_plot("Figure4_GS_MM_plots.pdf", pg, base_height = 6, base_width = 9)


nw1_net <- network.initialize(length(nw1_idx), directed = FALSE)
nw1_adj <- bicor(t(beta[nw1_idx,]), use = "p", maxPOutliers = 0.1)
nw1_adj[nw1_adj < 0.8] <- 0
nw1_adj <- nw1_adj^12
network.adjacency(nw1_adj, nw1_net, ignore.eval = FALSE, names.eval = "adj")
nw1$GeneNames[is.na(nw1$GeneNames)] <- "NA"
nw1_net %v% "GeneNames" <- nw1$GeneNames
cols_nw1 <- character(length(unique(nw1$GeneNames)))
names(cols_nw1) <- unique(nw1$GeneNames)
cols_nw1[1:8] <- c("grey", brewer.pal(7, "Paired"))
plot_nw1_net <- ggnet2(nw1_net, color = "GeneNames", edge.size = "adj", size = 3, palette = cols_nw1, legend.position = "right")



# ## RDA
# 
# library(vegan)
# library(ggplot2)
# library(MASS)
# library(qvalue)
# 
# rda_full <- rda(t(beta[idx,]) ~ ., phenos[,c(2,5)])
# rda_full_test <- anova.cca(rda_full)
# rda_full_R2 <- RsquareAdj(rda_full)
# 
# sc <- scores(rda_full, scaling = 3, display = c("sp", "wa", "bp"))
# 
# rda <- ggplot() +
#     theme_bw() +
#     theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 10, colour = "black"), legend.text = element_text(size = 10, colour = "black")) +
#     geom_hline(yintercept = 0, colour = "grey50", linetype = 3) +
#     geom_vline(xintercept = 0, colour = "grey50", linetype = 3) +
#     #geom_point(data = as.data.frame(scc$species), aes(x = RDA1, y = PC1), color = "grey32", size = 0.5, shape = "+") +
#     geom_point(data = as.data.frame(sc$sites), mapping = aes(x = RDA1, y = RDA2), size = 3) +
#     geom_text(data = as.data.frame(sc$sites), mapping = aes(x = RDA1, y = RDA2, label = rownames(as.data.frame(sc$sites))), nudge_x = 0.2, nudge_y = 0.2, size = 2) +
#     geom_segment(data = as.data.frame(sc$biplot), aes(x = 0, y = 0, xend = RDA1, yend = RDA2), size = 1, arrow = arrow(length = unit(0.03, "npc"))) +
#     geom_label(data = as.data.frame(sc$biplot), mapping = aes(x = RDA1, y = RDA2, label = row.names(sc$biplot)), fontface = "bold", size = 2, nudge_x = 0.2) +
#     labs(title = "A")
# rda
# 
# # center and scale scores
# resscale <- apply(as.matrix(sc$species[,1]), 2, scale)
# 
# qvals <- matrix(nrow = nrow(resscale), ncol = ncol(resscale), dimnames = list(rownames(resscale), "RDA1"))
# 
# # Calculate significance based on robust Mahalanobis distance for each variable independently
# for(i in 1:ncol(resscale)) {
#     rcov <- cov.rob(resscale[, i]) # robust covariance
#     resmaha <- mahalanobis(as.matrix(resscale[, i]), rcov$center, rcov$cov) # test statistic
#     lambda <- median(resmaha) / qchisq(0.5, df = 1)
#     pvals <- pchisq(resmaha / lambda, 1, lower.tail = FALSE)
#     qvals[, i] <- qvalue(pvals)$qvalues
# }
# 
# 
# # Manhattan plot of env_stp_rda score p-values
# qsig <- which(qvals < 0.05)
# regions <- getData(tiles)[wt_var > 0.05, 1:3]
# chrs <- regions$chr
# man <- ggplot(data.frame(pos = 1:length(pvals), pval = -log10(pvals)), aes(x = pos, y = pval)) +
#     geom_point(colour = c("grey", "grey30")[(as.numeric(chrs)%%2)+1], size = 0.1) +
#     ylab(expression("-log"[10]*"(pvalue)")) +
#     scale_x_continuous(name = "Chromosome", breaks = c(cumsum(table(chrs))-table(chrs)/2), labels =1:29) +
#     geom_point(data = data.frame(pos = (1:length(pvals))[qsig], pval = -log10(pvals)[qsig]), aes(x = pos, y = pval), color = "blue", size = 0.5) +
#     theme_bw() +
#     theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(size = 6), legend.text = element_text(size = 10, colour = "black")) +
#     labs(title = "C")
# ggsave(filename = "RDA_manhattan.png", man, "png", width = 8, height = 4, units = "in", dpi = 300)
# 
# rda_regions <- GRanges(seqnames = chrs[qsig], IRanges(start = regions$start[qsig], end = regions$end[qsig]))
# 
# hits <- findOverlaps(rda_regions, dmrs)
# 
# df <- data.frame(width = dmrs@ranges@width, overlap = as.numeric(1:length(dmrs) %in% subjectHits(hits)))
# ggplot(df, aes(x = width, fill = factor(overlap))) +
#     geom_density(alpha = 0.3)
# 
# ## Plotting highlighting SNPs
# col_snps <- rep("grey", ncol(all_freq))
# col_snps[which(colnames(all_freq) %in% env_rda_outliers_C)] <- "blue"
# empty <- col_snps
# empty[grep("grey30",empty)] <- rgb(0,1,0, alpha=0) # transparent
# empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","blue")
# 
# snps <- ggplot() +
#     theme_bw() +
#     theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 10, colour = "black"), legend.text = element_text(size = 10, colour = "black")) +
#     geom_hline(yintercept = 0, colour = "grey50", linetype = 3) +
#     geom_vline(xintercept = 0, colour = "grey50", linetype = 3) +
#     geom_point(data = as.data.frame(scc$species), aes(x = RDA1, y = PC1), color = col_snps, fill = col_snps, size = c(0.1,0.5)[(col_snps == "blue")+1], shape = 21) +
#     labs(title = "B")
# #    geom_segment(data = as.data.frame(scc$biplot/10), aes(x = 0, y = 0, xend = RDA1, yend = PC1), size = 1, arrow = arrow(length = unit(0.03, "npc"))) +
# #    geom_label(data = as.data.frame(scc$biplot/10), mapping = aes(x = RDA1, y = PC1, label = "Temp./Elev."), fontface = "bold", nudge_x = -0.02)
# snps
# 
