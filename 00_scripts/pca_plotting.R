#!/usr/bin/env Rscript
#require(adegenet, quietly = T)
require(RColorBrewer, quietly = T)

args <- commandArgs(trailingOnly = T)
# system("/Applications/plink_mac/plink --bfile ~/Desktop/CAST_2016_50K_analysis/01_data/CAST_2016_50K_neutral --chr-set 29 --pca --out ~/Desktop/CAST_2016_50K_analysis/01_data/CAST_2016_50K_neutral")
# args <- "~/Desktop/CAST_2016_50K_analysis/01_data/CAST_2016_50K_neutral"

if(length(args) == 0) {stop("You must provide a file prefix")}
if(length(args) > 1) {warning("More than one argument detected! Using the first argument as the file prefix")}

f <- paste0(args[1], ".eigenvec")
e <- paste0(args[1], ".eigenval")

if(!file.exists(f)) {stop("The file prefix you supplied was not found. Have you conducted the PCA yet?")}

d <- read.table(f, header = F)
d[,1] <- droplevels(as.factor(unlist(as.data.frame(strsplit(as.character(d[,1]), "_"))[1,])))
eig <- read.table(e, header = F)
eig <- round(eig[,1] / sum(eig[,1]) * 100, digits = 1)
mycol <- paste0(c(brewer.pal(12,"Paired"), "#000000", "#808080", "#CCCCCC"), "B3") #B3 = 70%


pdf(paste0(args[1], "_plot_axes1_2.pdf"), width = 3.5, height = 3.5)
par(mfrow = c(2,2), mar = c(0,0,0,0), oma = c(5,5,1,1), xpd = NA)
plot(x = d[,3], y = d[,4], 
     col = mycol[d[,1]],
     pch = 16,
     cex = 0.7,
     cex.axis = 0.7,
     ylab = paste0("PC2 (", eig[2], "%)"), 
     xlab = "",
     xaxt = "n")
plot(x = d[,5], y = d[,4], 
     col = mycol[d[,1]],
     pch = 16,
     cex = 0.7,
     cex.axis = 0.7,
     ylab = "", 
     yaxt = "n",
     xlab = paste0("PC3 (", eig[3], "%)"))
plot(x = d[,3], y = d[,5], 
     col = mycol[d[,1]],
     pch = 16,
     cex = 0.7,
     cex.axis = 0.7,
     ylab = paste0("PC3 (", eig[3], "%)"), 
     xlab = paste0("PC1 (", eig[1], "%)"))
plot.new()
legend("bottom", legend = levels(d[,1]), pch = 16, col = "grey", pt.bg = mycol, pt.cex = 0.5, cex = 0.5, ncol = 3)
dev.off()

#### manual plotting code

# anno <- read.table("~/Desktop/CAST_2016_50K_analysis/01_data/CAST_2016_50K_lowlink.fam")
# rownames(anno) <- as.character(anno$V2)
# anno <- anno[as.character(d$V2),]
# col <- funky(nlevels(d[,1]))[d[,1]] # by sub-basin
# col <- funky(2)[anno[,5]] # by sex
# col <- funky(2)[anno[,6]] # by CR good or bad
# plot(x = d[,3], y = d[,4], col = col, pch = 16, ylab = paste0("PC2 (", eig[2], "%)"), xlab = paste0("PC1 (", eig[1], "%)"))
# pch <- unlist(lapply(levels(d[,1]), function(i) {as.numeric(factor(gsub("^[A-Z].*_", "", gsub("_s[0-9].*", "", d[d[,1] == i, 2]))))}))
# plot(x = d[,3], y = d[,4], col = col, pch = pch, ylab = paste0("PC2 (", eig[2], "%)"), xlab = paste0("PC1 (", eig[1], "%)"))

# data <- read.PLINK("~/Desktop/CAST_2016_50K_analysis/01_data/CAST_2016_50K_lowlink.raw", quiet = T)
# mis <- sapply(NA.posi(data), length)/nLoc(data)
# col <- colorRampPalette(c("firebrick3", "yellow", "steelblue3"))(8)[cut(mis, seq(0, 0.08, 0.01))]
# plot(x = d[,3], y = d[,4], col = col, pch = 16, ylab = paste0("PC2 (", eig[2], "%)"), xlab = paste0("PC1 (", eig[1], "%)")); legend("topleft", levels(cut(mis, seq(0, 0.08, 0.01))), title = "Prop. missing data", pch = 16, col = colorRampPalette(c("firebrick3", "yellow", "steelblue3"))(8))

# selection of the three groups present in the full dataset
# abline(a=-0.04,b=-1.5, lty = 2)
# abline(a=0.06,b=-1.5, lty = 2)

# grp_vec <- numeric(nrow(d))
# grp_vec[d[,4] < ((-1.5*d[,3]) - 0.04)] <- 1
# grp_vec[d[,4] > ((-1.5*d[,3]) + 0.06)] <- 3
# grp_vec[grp_vec == 0] <- 2
# plot(x = d[,3], y = d[,4], col = grp_vec, pch = 16, ylab = paste0("PC2 (", eig[2], "%)"), xlab = paste0("PC1 (", eig[1], "%)"))
# write.table(cbind(d[,1:2],paste0("pop",grp_vec)), "~/Desktop/CAST_2016_50K_analysis/01_data/CAST_2016_50K_full_pca.clusters", row.names = F, col.names = F, quote = F)
# plot(hclust(dist(d[,3:4])))
