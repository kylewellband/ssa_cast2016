#!/usr/bin/env Rscript
# popfreq2bayescan.R
# by: Kyle Wellband
# Function takes a population MAF file as produced by PLINK using --freq and 
# either --family or --within <groupfile> options and converts to bayescan input

args <- commandArgs(trailingOnly = T)

freq.name <- args[1]
outfile.name <- sub(".frq.strat", ".bayescan", freq.name)

freq <- read.table(freq.name, header = T)

nloci <- nlevels(freq$SNP)

pops <- levels(freq$CLST)

bayescan_data <- list()

for(i in 1:length(pops)) {
  bayescan_data[[i]] <- matrix(nrow = nloci, ncol = 5)
  bayescan_data[[i]][,1] <- 1:nloci
  bayescan_data[[i]][,2] <- freq$NCHROBS[freq$CLST == pops[i]]
  bayescan_data[[i]][,3] <- 2
  bayescan_data[[i]][,4] <- freq$NCHROBS[freq$CLST == pops[i]] - freq$MAC[freq$CLST == pops[i]]
  bayescan_data[[i]][,5] <- freq$MAC[freq$CLST == pops[i]]
}

out_vec <- list()  # character(4+(3*length(pops))+(nloci*length(pops)))
  
out_vec[1:4] <-  c(paste0("[loci]=",nloci), "", paste0("[populations]=",length(pops)), "")

out_vec[seq(1,length(pops)*3, 3)+4] <- paste0("[pop]=", 1:length(pops))

out_vec[seq(1,length(pops)*3, 3)+6] <- ""

sub <- seq(1,length(pops)*3, 3)+5

for(i in 1:length(pops)) {
  out_vec[[sub[i]]] <- apply(bayescan_data[[i]], 1, paste, collapse = " ")
}

write.table(unlist(out_vec), outfile.name, row.names = F, col.names = F, quote = F)
