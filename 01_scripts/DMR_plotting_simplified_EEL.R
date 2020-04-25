library("data.table")
library("ComplexHeatmap")

ldmrs <- fread("DSS_results_pop_cov10_pop_DMR_heatmap_mean_Betas.txt")
sample_info <- read.table("sample_info_pop_indvcov10.txt", header = TRUE, stringsAsFactors = FALSE)
anno_columns <- "pop"
coef <- "pop"
design <- data.frame(sample_info[, anno_columns])
design[] <- lapply(design, factor)
names(design) <- anno_columns

col_cols <- lapply(anno_columns, function(i) {
    levels <- unique(sample_info[, i])
    coef_cols <- sapply(1:length(levels), function(j) {grey.colors(length(levels))[j]})
    names(coef_cols) <- levels
    coef_cols
})
names(col_cols) <- anno_columns

if (grepl(":", coef)) {
    col_order <- order(interaction(sample_info[ , anno_columns]))
    coef <- gsub(":", "\\.", coef)
} else {
    col_order <- order(interaction(sample_info[ , c(anno_columns[!anno_columns %in% coef], coef)]))
}

col_anno <- HeatmapAnnotation(df = design[col_order, ], col = col_cols)

png(filename = paste0("DSS_results_pop_cov10_pop_DMR_heatmap.png"), width = 8, height = 11, units = "in", res = 300)

print(Heatmap(
    matrix = as.matrix(ldmrs[,-1:-5])[, col_order],
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    row_title = NULL,
    cluster_columns = FALSE,
    use_raster = TRUE,
    raster_device = "png",
    top_annotation = col_anno
))

dev.off()

