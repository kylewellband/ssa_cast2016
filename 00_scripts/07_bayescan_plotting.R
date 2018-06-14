setwd("~/Desktop/CAST_2016_50K_analysis/")

source("/Applications/BayeScan2.1/R functions/plot_R.r")

map <- read.table("./01_data/CAST_2016_50K.map", stringsAsFactors = F)
remove.outliers <- plot_bayescan("./02_results/07_outliers/bayescan/pr_odds_10_fst.txt", add_text = F)$outliers
write.table(map[remove.outliers, 2], "./01_data/CAST_2016_50K_outliers_to_remove_for_neutral.txt", quote = F, row.names = F, col.names = F)

pdf("./02_results/07_outliers/bayescan/CAST_2016_50K_bayescan_pr_odd_10000.pdf", width = 8, height = 8)
bayescan.outliers <- plot_bayescan("./02_results/07_outliers/bayescan/pr_odds_10000_fst.txt", add_text = F)$outliers
dev.off()

write.table(map[bayescan.outliers, 2], "./02_results/07_outliers/bayescan_outliers.accnos", quote = F, row.names = F, col.names = F)

