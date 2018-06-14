setwd("~/Desktop/CAST_2016_50K_analysis/")

library(pcadapt)
library(qvalue)
library(qqman)

data <- read.pcadapt("01_data/CAST_2016_50K.ped", type = "ped")
annot <- read.table("01_data/CAST_2016_50K.fam")
snp_annot <- read.table("01_data/CAST_2016_50K.bim")

pca <- pcadapt(data, K = 2, min.maf = 0.01)

plot(pca, option = "screeplot")

plot(pca, option = "scores")

if (!dir.exists("./02_results/99_fusion/")) {
    dir.create("./02_results/99_fusion/")
}

df <- data.frame(CHR = snp_annot[, 1], POS = snp_annot[, 4], SNP = snp_annot[, 2], pvalues = pca$pvalues)
df$pvalues[which(df$pvalues==0)] <- min(df$pvalues[-which(df$pvalues==0)])

write.table(df, "./02_results/99_fusion/pcadapt_pvals.txt", quote = F)

pdf("./02_results/99_fusion/pcadapt_manhattan_plot.pdf", width = 3.5, height = 2, pointsize = 4)
#png("./02_results/99_fusion/pcadapt_manhattan_plot.png", width = 3.5, height = 2, units = "in", pointsize = 4, res = 300)
#manhattan(df, bp = "POS", p = "pvalues", logp = T, suggestiveline = F, genomewideline = F, ylab = expression("-log"[10]*"(pvalues)"), ylim = c(0,350), cex = 0.5, cex.axis = 1)
#dev.off()
ggplot(df, aes(x = 1:nrow(df), y = -log10(pvalues))) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"), legend.text = element_text(size = 12, colour = "black")) +
    geom_point(color = c("grey30", "grey")[(as.numeric(df$CHR)%%2) + 1], size = 0.01) +
    ylab(expression("-log"[10]*"(pvalues)")) +
    xlab("Chromosome") +
    scale_x_continuous(breaks = as.numeric(round(cumsum(table(df$CHR))-table(df$CHR)/2)), labels = c(rep("",7),8,rep("",20), 29))
dev.off()


library(patchwork)
