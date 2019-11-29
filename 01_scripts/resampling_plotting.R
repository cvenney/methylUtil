library("DSS")
library("data.table")
library("tidyverse")

dss_res <- fread("Projects/sasa_epi/03_results/DSS_DML_results.txt.gz")
DML <- nrow(callDML(dss_res, delta = 0.2, p.threshold = 0.01))
DMR <- nrow(callDMR(dss_res, delta = 0.2, p.threshold = 0.01))
rm(dss_res)

x<-fread("Projects/sasa_epi/03_results/resample_results.txt.gz")

loci <- rowSums(pivot_wider(x, id_cols = rep, names_from = chr, values_from = nloci)[,-1])

regions <- rowSums(pivot_wider(x, id_cols = rep, names_from = chr, values_from = nregion)[,-1])

ggplot(data.frame(DML = loci), aes(x = DML)) +
    geom_histogram(binwidth = 250) +
    geom_vline(xintercept = DML) +
    ylab("Count")

ggplot(data.frame(DMR = regions), aes(x = DMR)) +
    geom_histogram(binwidth = 25) +
    geom_vline(xintercept = DMR) +
    ylab("Count")
