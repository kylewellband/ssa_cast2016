#!/usr/bin/env Rscript

setwd("~/Desktop/CAST_2016_50K_analysis/")
sum_stat_prefix <- "./02_Results/00_Basic_Stats/"
het <- system(paste0("ls ", sum_stat_prefix, "*.hwe"), intern = T)
ho <- he <- numeric(length(het))
for(p in 1:length(het)) {
  d = read.table(het[p], header = F)
  ho[p] = mean(d[,7])
  he[p] = mean(d[,8])
  rm(d)
}

ibc <- system(paste0("ls ", sum_stat_prefix, "*.ibc"), intern = T)
ib <- numeric(length(het))
n <- numeric(length(het))
for(p in 1:length(het)) {
  d = read.table(ibc[p], header = T)
  ib[p] = mean(d[,4])
  n[p] = nrow(d)
  rm(d)
}
sum_stat_df <- data.frame(Pop = gsub(".hwe", "", gsub(sum_stat_prefix, "", het)),
                          N = n,
                          Ho = round(ho, digits = 3),
                          He = round(he, digits = 3),
                          Fis = round(ib, digits = 3))

write.table(sum_stat_df, "./02_Results/00_Basic_Stats/00_sum_stat_table.txt", quote = F, row.names = F, sep = "\t")
