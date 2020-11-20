#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

prob <- 0.9

bg <- fread(args[1], header = FALSE)

bg <- bg[V1 %like% "NC_0[1-9].*", ]

bg[, cov := V5 + V6]

bg <- bg[cov <= quantile(bg$cov, 0.995),]

pval <- -log10(pbinom(q = bg$V5, size = bg$cov, prob = prob))

out <- cbind(bg[, 1:3], pval, bg[,cov])

fwrite(out, paste0("05_bed_files/combined_coverage_logpval_lambda", prob, ".bedGraph.gz"), sep = "\t", col.names = FALSE)
