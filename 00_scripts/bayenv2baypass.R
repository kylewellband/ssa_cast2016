#!/usr/bin/env Rscript

setwd("~/Desktop/CAST_2016_50K_analysis/")

x <- as.matrix(read.table("./01_data/CAST_2016_50K.bayenv"))

z <- matrix(ncol = ncol(x)*2, nrow = nrow(x)/2)

index <- 1
ref_allele <- seq(1, ncol(x)*2, 2)
alt_allele <- seq(2, ncol(x)*2, 2)

allele <- 0

for (i in 1:nrow(x)) {
    if (allele == 0) {
        z[index, ref_allele] <- x[i,]
        allele <- 1
    } else {
        z[index, alt_allele] <- x[i,]
        allele <- 0
        index <- index + 1
    }
    
}

write.table(z, "./01_data/CAST_2016_50K.baypass", col.names = F, row.names = F, quote = F)
