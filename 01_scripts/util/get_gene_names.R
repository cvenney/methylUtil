#!/usr/bin/env Rscript 
# get_gene_names.R

args <- commandArgs(TRUE)
# setwd("~/Projects/sasa_epi/methylUtil"); args <- '06_methylation_results/adults_6x6_min5_max20_groupWild_dmr_pval0.001_geneids.txt'

if (length(args) != 1) {
    stop("USAGE: get_gene_name.R <gene ID file>")
}

if (!file.exists("02_reference/feature_table.txt.gz")) {
    stop("File: 02_reference/feature_table.txt.gz, does not exist!")
}

if (!suppressPackageStartupMessages(require(data.table))) {
    install.packages("data.table")
    suppressPackageStartupMessages(library(data.table))
}

ref <- fread("02_reference/feature_table.txt.gz")
ref <- ref[feature == "mRNA"]
ref[, LOC := as.character(GeneID)]

geneid <- fread(args[1], header = FALSE)
if (any(grepl("^gene", geneid$V1))) {
    geneid[, LOC := sub("gene[0-9]+:", "", V1)]
} else {
    geneid[, LOC := V1]
}

geneid <- merge(x = geneid, y = ref, by.x = "LOC", by.y = "LOC", all.x = TRUE)
geneid[, NAME := sub(", transcript variant X[0-9]+", "", name)]

names <- unique(geneid[,.("GeneID" = V1, "Symbol" = symbol, NAME)])

fwrite(names, file = sub("geneids.txt", "gene_names.txt", args[1]), sep = "\t")