# FLK outlier analysis

setwd("~/Desktop/CAST_2016_50K_analysis/")
source("./00_scripts/FLK.R")
library(reshape2)
library(ape)

# Create population allele frequency file for neutral and all SNPs
system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K_neutral --chr-set 29 --freq --family --out ./01_data/CAST_2016_50K_neutral")
system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --freq --family --out ./01_data/CAST_2016_50K")

# read allele frequency data
ll_frq <- read.table("./01_data/CAST_2016_50K_neutral.frq.strat", header = T) # neutral SNPs
all_frq <- read.table("./01_data/CAST_2016_50K.frq.strat", header = T) # all SNPs

# set up output directory
dir.create("./02_results/07_outliers/FLK")
setwd("./02_results/07_outliers/FLK")

# reorganize allele frequency data into matrices of pops (rows) X SNP freqs (columns)
ll_freq <- dcast(ll_frq, CLST ~ SNP, value.var = "MAF")
rownames(ll.freq) <- ll_freq[,1]
ll_freq <- ll_freq[,-1]

all_freq <- dcast(all_frq, CLST ~ SNP, value.var = "MAF")
rownames(all_freq) <- all_freq[,1]
all_freq <- all_freq[,-1]

# Calculate demographic relationship among pops using neutral SNPS
Ftree <- Fij(ll_freq, outgroup = "NWM") # this take a couple hours but is saved as fij.txt

# Conduct FLK tests based on Chi-Sq distribution
tests <- FLK(all_freq, Ftree)

# Adjust for multiple testing using FDR (could alternatively use qvalue package)
sum(p.adjust(tests$F.LK.p.val)<0.05)
rownames(tests)[which(p.adjust(tests$F.LK.p.val)<0.05)]
write.table(rownames(tests)[which(p.adjust(tests$F.LK.p.val)<0.05)], "../FLK_outliers.accnos", col.names = F, row.names = F, quote = F)

# run simulations to estimate the neutral envelope
#system("python2 ../../../00_scripts/FLKnull.py 10000")

# make FLK plot
source("../../../00_scripts/plotNull.R")

points(tests$Ht, tests$F.LK, pch = 16, cex = 0.3)

col <- rep(1,nrow(tests));col[which(p.adjust(tests$F.LK.p.val)<0.05)] <- 2
pch <- rep(16,nrow(tests));pch[which(p.adjust(tests$F.LK.p.val)<0.05)] <- 8
pcex <- rep(0.3,nrow(tests));pcex[which(p.adjust(tests$F.LK.p.val)<0.05)] <- 1
plot(tests$Ht, tests$F.LK, col = col, pch = pch, cex = pcex)

# save all results
write.table(tests, "../FLK_results.txt", sep = "\t", quote = F)
