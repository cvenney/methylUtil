#!/usr/bin/env Rscript

#setwd("~/Projects/sasa_epi/methylUtil/")
## Install and load necessary packages
for (p in c("data.table", "tidyverse", "circlize", "ComplexHeatmap", "ggplot2", "ggforce", "GenomicRanges")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if (p %in% c("GenomicRanges")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

col_cols <- c("gold", "blue")
names(col_cols) <- c("SAS", "Wild")

age_cols <- c("black", "grey")
names(age_cols) <- c("Adult", "Juvenile")

lg1 <- Legend(col_fun = colorRamp2(c(0, 1), c("white", "red")),
              title = "Methylation", border = "black")
lg2 <- Legend(labels = names(col_cols), legend_gp = gpar(fill = col_cols), title = "Origin")
lg3 <- Legend(labels = names(age_cols), legend_gp = gpar(fill = age_cols), title = "Age")
pd <- packLegend(lg1, lg2, lg3)

adult <- fread("06_methylation_results/adults_6x6_group_DMR_heatmap_mean_Betas.txt")
asamples <- data.table(sampleID = names(adult)[-1:-5])
asamples[, Group := ifelse(grepl("W", sampleID), "Wild", "SAS")]
asamples[, Age := "Adult"]
acol_anno <- HeatmapAnnotation(df = asamples[, .(Group, Age)], 
                               col = list(Group = col_cols, Age = age_cols), 
                               show_annotation_name = FALSE,
                               simple_anno_size = unit(3, "mm"))

a_hp <- Heatmap(
    matrix = as.matrix(adult[,-1:-5]),
    name = "Methylation", 
    border = "black",
    col = colorRamp2(c(0, 1), c("white", "red")),
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    row_title = NULL,
    cluster_columns = TRUE,
    clustering_distance_columns = "pearson",
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_device = "png",
    top_annotation = acol_anno
)


juvenile <- fread("06_methylation_results/juveniles_8x8_group_DMR_heatmap_mean_Betas.txt")
jsamples <- data.table(sampleID = names(juvenile)[-1:-5])
jsamples[, Group := ifelse(grepl("W", sampleID), "Wild", "SAS")]
jsamples[, Age := "Juvenile"]
jcol_anno <- HeatmapAnnotation(df = jsamples[, .(Group, Age)], 
                               col = list(Group = col_cols, Age = age_cols), 
                               annotation_legend_param = list(direction = "vertical"), 
                               show_annotation_name = FALSE, show_legend = FALSE,
                               simple_anno_size = unit(3, "mm"))
jdend <- as.dendrogram(hclust(as.dist(1 - cor(t(juvenile[, -1:-5]), method = "pearson")), method = "complete"))
j_hp <- Heatmap(
    matrix = as.matrix(juvenile[, -1:-5]),
    name = "Methylation",
    border = "black",
    col = colorRamp2(c(0, 1), c("white", "red")),
    cluster_rows = jdend,
    row_title = NULL,
    cluster_columns = TRUE,
    clustering_distance_columns = "pearson",
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_device = "png",
    top_annotation = jcol_anno,
    heatmap_legend_param = list(direction = "vertical", border = "black"),
    show_heatmap_legend = FALSE
)

## Common regions heatmap
aGR <- GRanges(seqnames = adult$seqnames, ranges = IRanges(start = adult$start, end = adult$end))
jGR <- GRanges(seqnames = juvenile$seqnames, ranges = IRanges(start = juvenile$start, end = juvenile$end))

hits <- findOverlaps(aGR, jGR)

csamples <- rbind(asamples, jsamples)
ccol_anno <- HeatmapAnnotation(df = csamples[order(csamples$Group, csamples$Age), .(Group, Age)], 
                               col = list(Group = col_cols, Age = age_cols), 
                               annotation_legend_param = list(direction = "vertical"), 
                               show_annotation_name = FALSE,
                               simple_anno_size = unit(3, "mm"))
cmat <- as.matrix(cbind(adult[queryHits(hits),-1:-5], juvenile[subjectHits(hits),-1:-5]))[, order(csamples$Group, csamples$Age)]
cmat[3,] <- colMeans(cmat[3:5,])
c_hp <- Heatmap(
    matrix = cmat,
    name = "Methylation",
    border = "black",
    col = colorRamp2(c(0, 1), c("white", "red")),
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    row_title = NULL,
    cluster_columns = TRUE, 
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_device = "png",
    top_annotation = ccol_anno,
    heatmap_legend_param = list(direction = "vertical", border = "black"),
    show_heatmap_legend = FALSE
)


# ## Venn of overlap
# df_vdc <- data.frame(Counts = c("Adult\n260", "Juvenile\n359", "3")) %>%
#     mutate(x = c(-1, 1, 0), y = c(0, 0, 0))
# 
# df_venn <- data.frame(x = c(-0.5, 0.5),
#                       y = c(0, 0),
#                       labels = c("Adult", "Juvenile"))
# 
# venn_plot <- ggplot(df_venn) +
#     geom_circle(aes(x0 = x, y0 = y, r = 1, fill = labels), alpha = 0.5, size = 1, colour = 'black', show.legend = FALSE) +
#     coord_fixed() +
#     theme_void() +
#     labs(fill = NULL) +
#     annotate("text", x = df_vdc$x, y = df_vdc$y, label = df_vdc$Counts, size = 4) +
#     scale_fill_brewer(palette = "Dark2")

## Set up viewports

adult_vp <- viewport(x = unit(0.25, "npc"), y = unit(0.66, "npc"),
                     width = unit(0.4, "npc"), height = unit(0.66, "npc"))

juvenile_vp <- viewport(x = unit(0.75, "npc"), y = unit(0.5, "npc"),
                        width = unit(0.5, "npc"), height = unit(1, "npc"))

# venn <- viewport(x = unit(0.25, "npc"), y = unit(0.17, "npc"),
#                  width = unit(0.5, "npc"), height = unit(0.25, "npc"))

common_vp <- viewport(x = unit(0.25, "npc"), y = unit(0.17, "npc"),
                 width = unit(0.4, "npc"), height = unit(0.25, "npc"))

## Actually draw the figure

pdf("Final_heatmap.pdf", width = 8, height = 8)
grid.newpage()

pushViewport(adult_vp)
draw(a_hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, newpage = FALSE)
popViewport()

pushViewport(viewport(x = unit(0.03, "npc"), y = unit(0.95, "npc"), width = unit(0.1, "npc"), height = unit(0.1, "npc")))
grid.text("A", gp = gpar(cex = 2, fontface = "bold"))
popViewport()

pushViewport(juvenile_vp)
draw(j_hp, annotation_legend_list = pd, newpage = FALSE)
popViewport()

pushViewport(viewport(x = unit(0.5, "npc"), y = unit(0.95, "npc"), width = unit(0.1, "npc"), height = unit(0.1, "npc")))
grid.text("B", gp = gpar(cex = 2, fontface = "bold"))
popViewport()

# pushViewport(venn)
# print(venn_plot, newpage = FALSE)
# popViewport()

pushViewport(common_vp)
draw(c_hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, newpage = FALSE)
popViewport()

pushViewport(viewport(x = unit(0.03, "npc"), y = unit(0.25, "npc"), width = unit(0.1, "npc"), height = unit(0.1, "npc")))
grid.text("C", gp = gpar(cex = 2, fontface = "bold"))
popViewport()

dev.off()
