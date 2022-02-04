#!/usr/bin/env Rscript
#srun -c 1 --mem 50G -p large --time 21-00:00:00 -J DSS_jackknife -o DSS_jackknife_%j.log Rscript ./01_scripts/07_brook_charr_DSS_jackknife.R &

library(DSS)
library(data.table)

BSobj <- readRDS(file = "06_methylation_results/brook_charr_lab_sample_info_all_data.rds")
files <- list.files(path = "samples_jackknife_family", pattern = ".txt")

for (i in 1 : length(files))
{

# samples and experimental design
print(files[i])
sample_info <- read.table(paste0("samples_jackknife_family/",files[i]), header=TRUE)
vars <- sample_info[,1:5]
#maternal jackknifing
#lab <- strsplit(files[i], c("_"))[[1]][1]
#family jackknifing
lab <- strsplit(files[i], c("\\."))[[1]][1]


#subset <- BSobj[1:1000, sample_info$sample]
subset <- BSobj[, sample_info$sample]					#remove 1:1000* later - for subset

### DML testing
DMLfit = DMLfit.multiFactor(subset, design=vars, formula=~adult_temp + juv_temp + adult_temp:juv_temp)
colnames(DMLfit$X)

# adult temp
DMLtest.ad = DMLtest.multiFactor(DMLfit, coef="adult_tempfroid")

	ix=sort(DMLtest.ad[,"pvals"], index.return=TRUE)$ix
	head(DMLtest.ad[ix,])	
write.table(DMLtest.ad[ix,], file=(paste0("07_jackknife_results/brook_charr_all_CpGs_", lab, "_adult_tempfroid_DMLs.txt")))

DMRtest = callDMR(DMLtest.ad, p.threshold=0.05)
write.table(DMRtest, file=(paste0("07_jackknife_results/brook_charr_all_CpGs_", lab, "_adult_tempfroid_DMRs_p0.05.txt")))


# juvenile temp
DMLtest.ad = DMLtest.multiFactor(DMLfit, coef="juv_tempfroid")

	ix=sort(DMLtest.ad[,"pvals"], index.return=TRUE)$ix
	head(DMLtest.ad[ix,])	
write.table(DMLtest.ad[ix,], file=(paste0("07_jackknife_results/brook_charr_all_CpGs_", lab, "_juv_tempfroid_DMLs.txt")))

DMRtest = callDMR(DMLtest.ad, p.threshold=0.05)
write.table(DMRtest, file=(paste0("07_jackknife_results/brook_charr_all_CpGs_", lab, "_juv_tempfroid_DMRs_p0.05.txt")))


# adult x juvenile temp
DMLtest.ad = DMLtest.multiFactor(DMLfit, coef="adult_tempfroid:juv_tempfroid")

	ix=sort(DMLtest.ad[,"pvals"], index.return=TRUE)$ix
	head(DMLtest.ad[ix,])	
write.table(DMLtest.ad[ix,], file=(paste0("07_jackknife_results/brook_charr_all_CpGs_", lab, "_int_adult_tempfroid_juv_tempfroid_DMLs.txt")))

DMRtest = callDMR(DMLtest.ad, p.threshold=0.05)
write.table(DMRtest, file=(paste0("07_jackknife_results/brook_charr_all_CpGs_", lab, "_int_adult_tempfroid_juv_tempfroid_DMRs_p0.05.txt")))



}