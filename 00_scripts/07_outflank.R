#### OutFLANK outlier analyses ####

setwd("~/Desktop/CAST_2016_50K_analysis/")

library(OutFLANK)

system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --recode A --out ./01_data/CAST_2016_50K")

dat <- read.table("./01_data/CAST_2016_50K.raw", header = T, stringsAsFactors = F)

geno <- as.matrix(dat[, -1:-6])
geno[is.na(geno)] <- 9

fst.mat <- MakeDiploidFSTMat(SNPmat = geno, locusNames = gsub(pattern = "_[0-9]", "", names(dat)[-1:-6]), popNames = dat[,1])

res <- OutFLANK(FstDataFrame = fst.mat, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin = 0.1, NumberOfSamples = 15, qthreshold = 0.05)

pdf("./02_results/07_outliers/OutFLANK_distn.pdf", height = 4, width = 4)
OutFLANKResultsPlotter(res, withOutliers = T, NoCorr = T, Hmin = 0.1, binwidth = 0.005, Zoom = F, RightZoomFraction = 0.05, titletext = NULL)
dev.off()

write.table(res$results, "./02_results/07_outliers/OutFLANK.res", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(res$LocusName[res$OutlierFlag], "./02_results/07_outliers/OutFLANK_outliers.accnos", col.names = F, row.names = F, quote = F)