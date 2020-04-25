#!/usr/bin/env Rscript
## get_gene_names.R

args <- commandArgs(TRUE)
# setwd("~/Desktop/sasa_epi/methylUtil/"); args <- "06_methylation_results/juveniles_8x8_min5_max20_groupWild_dmr_pval0.001_dmr_context.txt"

if (length(args) != 1) {
    stop("Usage: get_gene_names.R <bedtools_intersect_outfile>")
}

if (suppressPackageStartupMessages(!require(data.table))) {
    install.packages("data.table");suppressPackageStartupMessages(require(data.table))
}

ft <- fread("02_reference/feature_table.txt.gz", sep = "\t")

x <- fread(args[1], sep = "\t")
x <- x[V9 != "region"]
x[, geneid := as.integer(sub("(;|,).*", "", sub(".*GeneID:", "", V15)))]
x <- unique(x[, .(V1, V2, V3, V4, V5, geneid)])

out <- merge(x, ft[feature != "gene"], by.x = "geneid", by.y = "GeneID", all.x = TRUE)
out <- unique(out[, .(V1, V2, V3, V4, V5, geneid, name)])

fwrite(out, paste0(sub("(.*)?\\..*$", "\\1", args[1]), "_gene_names.txt"), quote = FALSE, sep = "\t")


# #### Compare genes detected in Ssal
# europe <- fread("02_reference/Rodriguez-Barreto_genes.txt")
# europe <- unique(merge(europe, ft, by.x = "geneids", by.y = "symbol")[,. (geneids, GeneID)])
# any(europe$GeneID %in% x$geneid)
# ft[symbol == europe[GeneID %in% x$geneid]$geneids]
