#### Outlier comparisons ####

setwd("~/Desktop/CAST_2016_50K_analysis/")

bayescan.outliers <- read.table("./02_results/07_outliers/bayescan_outliers.accnos", stringsAsFactors = F)[,1]
outflank.outliers <- read.table("./02_results/07_outliers/OutFLANK_outliers.accnos", stringsAsFactors = F)[,1]
baypass.ranked <- read.table("./02_results/07_outliers/baypass_XtX_order.accnos", stringsAsFactors = F)[,1]
baypass.outliers <- read.table("./02_results/07_outliers/baypass_outliers.accnos", stringsAsFactors = F)[,1]
FLK.outliers <- read.table("./02_results/07_outliers/FLK_outliers.accnos", stringsAsFactors = F)[,1]

which(baypass.ranked %in% bayescan.outliers)
which(baypass.ranked %in% outflank.outliers)
which(baypass.ranked %in% FLK.outliers)

common.outliers <- intersect(intersect(bayescan.outliers, outflank.outliers), FLK.outliers)

length(common.outliers)

which(baypass.ranked %in% common.outliers)

map <- read.table("./01_data/CAST_2016_50K.bim", stringsAsFactors =F)

bed.sub <- which(map[,2] %in% common.outliers)

common.outliers.bed <- data.frame(chr = paste0("chr", map[bed.sub, 1]),
                                  start = map[bed.sub, 4] - 1,
                                  end = map[bed.sub, 4],
                                  name = map[bed.sub, 2])
write.table(common.outliers.bed, "./02_results/07_outliers/common_outliers.bed", sep = "\t", col.names = F, row.names = F, quote = F)


env.rda <- read.table("./02_results/07_outliers/ENV_RDA_SNP_scores.txt", header = T, stringsAsFactors = F)
env.outliers <- row.names(env.rda)[env.rda$qval < 0.05]

sum(env.outliers %in% bayescan.outliers)
sum(env.outliers %in% outflank.outliers)
sum(env.outliers %in% FLK.outliers)
sum(env.outliers %in% common.outliers)


which(baypass.ranked %in% env.outliers)

# Manhattan plot
library(RColorBrewer)
fst <- read.table("./01_data/CAST_2016_50K.fst", header = T, colClasses = c("numeric", "character", "numeric","numeric","numeric"))
axis_lab <- (table(fst$CHR)/2) + c(0, cumsum(table(fst$CHR))[-29])

pdf("./02_results/FigureSX_manhattan_plot.pdf", width = 8, height = 4)
plot(fst$FST, pch = 19, cex = 0.3, col = c("grey50", "black")[(as.numeric(fst$CHR) %% 2)+1], xaxt = "n",
     ylab = expression("F"[ST]), xlab = "")
axis(1, c(0,cumsum(table(fst$CHR))), labels = F)
axis(1, axis_lab, tick = F, labels = 1:29, las = 3, cex.axis = 0.5)
points(which(as.character(fst$SNP) %in% bayescan.outliers), fst$FST[as.character(fst$SNP) %in% bayescan.outliers], col = "red")
points(which(as.character(fst$SNP) %in% outflank.outliers), fst$FST[as.character(fst$SNP) %in% outflank.outliers], col = "blue", pch = 2)
points(which(as.character(fst$SNP) %in% FLK.outliers), fst$FST[as.character(fst$SNP) %in% FLK.outliers], col = "green", pch = 0)
legend("topright", pch = c(1,2,0), col = c("red", "blue", "green"), legend = c("Bayescan", "Outflank", "FLK"))
dev.off()
