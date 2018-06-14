# Build a plot for inversion results

setwd("~/Desktop/CAST_2016_50K_analysis/")

library(grid)
library(ggplot2)


#### Run PCA and plot axes 1 and 2 ####
system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --pca --out ./01_data/CAST_2016_50K")

d <- read.table("./01_data/CAST_2016_50K.eigenvec", header = F)
d[,1] <- droplevels(as.factor(unlist(as.data.frame(strsplit(as.character(d[,1]), "_"))[1,])))
eig <- read.table("./01_data/CAST_2016_50K.eigenval", header = F)
eig <- round(eig[,1] / sum(eig[,1]) * 100, digits = 1)
fus_gt <- read.table("./01_data/CAST_2016_50K_inversion_status.txt")
col <- c("#1b9e77","#d95f02","#7570b3")

#par(mar = c(5,4,1,1)+0.1)
#pdf("./02_results/99_fusion/CAST_2016_50K_all_snp_pca.pdf", width = 3, height = 3)
pca <- ggplot(d, aes(x = d[,3], y = d[,4])) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          plot.title = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 8, colour = "black"), 
          legend.text = element_text(size = 10, colour = "black")) +
    geom_point(colour = c("#1b9e77","#d95f02","#7570b3")[fus_gt$V3], size = 1) +
    ylab(paste0("PC2 (", eig[2], "%)")) +
    xlab(paste0("PC1 (", eig[1], "%)")) +
    labs(title = "A")
pca
#dev.off()

#### Manhattan plot ####

df <- read.table("./02_results/99_fusion/pcadapt_pvals.txt", header = T)

system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --within ./01_data/CAST_2016_50K_inversion_status.txt --fst --out ./02_results/99_fusion/karyotype")
df <- read.table("./02_results/99_fusion/karyotype.fst", header = T)

man <- ggplot(df, aes(x = 1:nrow(df), y = FST)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 8, colour = "black", margin = margin()), 
          axis.title.y = element_text(size = 8, colour = "black", margin = margin()),
          axis.title.x = element_blank()) +
    geom_point(color = c("grey", "grey30")[(as.numeric(df$CHR)%%2) + 1], size = 0.1) +
    ylab(expression("F"[ST])) +
    xlab("Chromosome") +
    scale_x_continuous(breaks = as.numeric(round(cumsum(table(df$CHR))-table(df$CHR)/2)), labels = c(rep("",7),"chr8",rep("",20), "chr29"))
man
ggsave(filename = "./02_results/99_fusion/pcadapt_manhattan_plot.png", man, device = "png", width = 60, height = 30, units = "mm", dpi = 300)


# inset size
man <- ggplot(df, aes(x = 1:nrow(df), y = -log10(pvalues))) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.ticks = element_line(size = 0.1), 
          axis.text = element_text(size = 3, colour = "black", margin = margin()), 
          axis.title.y = element_text(size = 3, colour = "black", margin = margin()),
          axis.title.x = element_blank(),
          axis.ticks.length = unit(0.2, "mm"),
          plot.margin = margin(1,1,1,1)) +
    geom_point(color = c("grey", "grey30")[(as.numeric(df$CHR)%%2) + 1], size = 0.001) +
    ylab(expression("-log"[10]*"(pvalues)")) +
    xlab("Chromosome") +
    scale_x_continuous(breaks = as.numeric(round(cumsum(table(df$CHR))-table(df$CHR)/2)), labels = c(rep("",7),"chr8",rep("",20), "chr29"))
man
ggsave(filename = "./02_results/99_fusion/pcadapt_manhattan_plot_inset.png", man, device = "png", width = 30, height = 15, units = "mm", dpi = 300)




#### Binning function

binning_function <- function(LDmat = NULL, fst = NULL, distances = NULL, binwidth = NULL, nbins = NULL) {
    if (is.null(LDmat)) stop("You must provide a matrix of pairwise LD")
    if (!is.matrix(LDmat) || ncol(LDmat) != nrow(LDmat)) stop("You must provide a square matrix of pairwise LD")
    res <- list()
    nsnps <- nrow(LDmat)
    if (!is.null(nbins) && is.null(binwidth)) {
        bin_names <- cut(1:nsnps, nbins)
        map_out <- tapply(distances, bin_names, max, na.rm = T)
    }
    if (!is.null(binwidth)) {
        bin_names <- ggplot2::cut_interval(distances, length = binwidth)
        nbins <- nlevels(bin_names)
        map_out <- tapply(distances, bin_names, max, na.rm = T)
    }
    res[["map"]] <- map_out
    if (!is.null(fst)) {
        fst_out <- tapply(fst, factor(bin_names), max, na.rm = T)
        res[["fst"]] <- fst_out
    }
    binned_mat <- matrix(nrow = nbins, ncol = nbins)
    for(i in 1:nlevels(bin_names)) {
        row_sub <- which(as.numeric(bin_names) == i)
        for (j in 1:nlevels(bin_names)) {
            col_sub <- which(as.numeric(bin_names) == j)
            binned_mat[i,j] <- ifelse(length(row_sub) == 0 || length(col_sub) == 0 , NA, max(LDmat[row_sub, col_sub], na.rm = T))
        }
    }
    res[["binnedLD"]] <- binned_mat
    return(res)
}

## END bin function

library(LDheatmap)
library(RColorBrewer)

map <- read.table("./01_data/chr8_29.bim")
n8 <- sum(map[, 1] == 8)
n29 <- sum(map[, 1] == 29)
map <- map[c(n8:1, (n8+1):(n8+n29)),]
len_chr8 <- 26425719
#map[1:n8, 4] <- len_chr8 - map[1:n8, 4]
map[(n8+1):nrow(map), 4] <- map[(n8+1):nrow(map), 4] + len_chr8
ld <- read.table("./01_data/chr8_29.ld")

#write.table(map[,c(2,4)], "./01_data/chr8_29_haploview.bim", col.names = F, row.names = F, quote = F)

LDmat <- as.matrix(ld)[c(n8:1, (n8+1):(n8+n29)), c(n8:1, (n8+1):(n8+n29))]

LDbin_100bins <- binning_function(LDmat, nbins = 100, distances = map[,4])

bin_width <- 1000000

LDbin <- binning_function(LDmat, binwidth = bin_width, distances = map[,4])
chr_split <- round(len_chr8/bin_width)
LDbin_with_break <- matrix(nrow = (nrow(LDbin$binnedLD)+1), ncol = (ncol(LDbin$binnedLD)+1))
LDbin_with_break[c(1:chr_split, (chr_split+2):nrow(LDbin_with_break)), c(1:chr_split, (chr_split+2):nrow(LDbin_with_break))] <- LDbin$binnedLD
LDmap_with_break <- numeric(length(LDbin$map)+1)
LDmap_with_break[c(1:chr_split, (chr_split+2):nrow(LDbin_with_break))] <- LDbin$map
LDmap_with_break[chr_split+1] <- LDmap_with_break[chr_split]

chr8_reorder <- c(chr_split:1, (chr_split+1):length(LDmap_with_break))

pdf("./02_results/99_fusion/LD_figure.pdf", width = 3.5, height = 3.5)
ll <- LDheatmap(LDbin_with_break[chr8_reorder, chr8_reorder], genetic.distances = LDmap_with_break[chr8_reorder], distances = "physical", LDmeasure = "r", title = "", flip = T, add.map = F)
dev.off()


#### Chromosome Models #### 
library(Gviz)

chr_mod <- data.frame(chrom = paste0("chr", c(8, 8, 8, 8, 29, 29)),
                 chromStart = c(0, 12000000, 13000000, 14000000, 0, 1000000),
                 chromEnd = c(12000000, 13000000, 14000000, 39000000, 1000000, 43400000),
                 name = c("rDNA","p8","q8.1","q8.2","q29.0","q29.1"),
                 gieStain = c("gpos75", "acen", "acen", "gneg", "acen", "gneg"))

test <- new("IdeogramTrack")
test@bandTable <- chr_mod

for(chr in c(8,29)) {
    chromosome(test) <- chr
    axis <- GenomeAxisTrack()
    len <- ifelse(chr == 8, yes = 26000000, no = 43400000)
    pdf(paste0("./02_results/99_fusion/chr", chr, "_model.pdf"), width = 3, height = 1)
    plotTracks(test, from = 0, to = 0, lwd = 0, col = "black")
    dev.off()
}



#### FST distribution ####
# Calculate FST for each marker
system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --fst --family --out ./01_data/CAST_2016_50K")
system("/Applications/plink_mac/plink --file ./01_data/fus_inv --chr-set 29 --fst --family --out ./01_data/fus_inv")

# Read FST
snps <- read.table("./01_data/CAST_2016_50K.fst", header = T)
fus_inv <- read.table("./01_data/fus_inv.fst", header = T)

# Plot distribution
pdf("./02_results/99_fusion/fusion_inversion_subbasin_divergence.pdf", width = 3, height = 2.5)
fst<-ggplot(data = snps, aes(x = snps$FST)) +
    geom_density(fill = "grey") +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 10, colour = "black", face = "plain"), legend.text = element_text(size = 10, colour = "black")) +
    ylim(c(0,100)) +
    geom_line(data = data.frame(y = c(0, 50), x = rep(fus_inv$FST[1], 2)), aes(x = x, y = y)) +
    #geom_line(data = data.frame(y = c(0, 75), x = rep(fus_inv$FST[2], 2)), aes(x = x, y = y)) +
    geom_point(data = data.frame(y = c(50), x = fus_inv$FST[1]), aes(x = x, y = y)) +
    geom_vline(xintercept = quantile(snps$FST, 0.95), lty = 2) + 
    xlab(expression("F"[ST])) +
    ylab("SNP Density") +
    labs(title = "D")
dev.off()    

#### diversity ####
library(tidyverse)
map <- read.table("./01_data/chr8_29.bim")
n8 <- sum(map[, 1] == 8)
n29 <- sum(map[, 1] == 29)
map <- map[c(n8:1, (n8+1):(n8+n29)),]
len_chr8 <- 26425719
map[1:n8, 4] <- len_chr8 - map[1:n8, 4]
map[(n8+1):nrow(map), 4] <- map[(n8+1):nrow(map), 4] + len_chr8 + 25000000
inversion_region <- data.frame(chr = c(8,8,29,29), x = c(19175644, len_chr8, len_chr8 + 25000000, 29484452 + 25000000), ymin = c(0.5,0.5,0.5,0.5), ymax = c(1,1,1,1))

max_chr8 <- max(map[map[,1] == 8, 4])
min_chr29 <- min(map[map[,1] == 29, 4])
scale <- floor(min_chr29 / max_chr8) - 1
data.sum$value[data.sum$mask == 1] = data.sum$value[data.sum$mask == 1] / scale
step <- 10
low.end <- max(data.sum$value[data.sum$name != "desert"])
up.start <- ceiling(max(data.sum$value[data.sum$name != "desert"]))
breaks <- seq(0, max(data.sum$value), step)
labels <- seq(0, low.end+step, step)
labels <- append(labels, scale * seq(from=ceiling((up.start + step) / step) * step, length.out=length(breaks) - length(labels), by=step))

col <- c("#1b9e77","#d95f02","#7570b3")

system("./00_scripts/99_karyotype_diversity.sh")

# Het
het.aa <- read.table("./01_data/chr8_29_AA.hwe", header = T)
het.aa <- het.aa[het.aa[,"TEST"] == "ALL",]
het.aa <- het.aa[c(n8:1, (n8+1):(n8+n29)),]
het.ab <- read.table("./01_data/chr8_29_AB.hwe", header = T)
het.ab <- het.ab[het.ab[,"TEST"] == "ALL",]
het.ab <- het.ab[c(n8:1, (n8+1):(n8+n29)),]
het.bb <- read.table("./01_data/chr8_29_BB.hwe", header = T)
het.bb <- het.bb[het.bb[,"TEST"] == "ALL",]
het.bb <- het.bb[c(n8:1, (n8+1):(n8+n29)),]

df_het <- data.frame(chr = map[,1],
                     pos = map[,4],
                     "AA" = het.aa$"O.HET.",
                     "AB" = het.ab$"O.HET.",
                     "BB" = het.bb$"O.HET.")
td_het <- gather(df_het, "Karyotype", "Ho", -chr, -pos, factor_key = T)
td_het$K2 <- factor(td_het$Karyotype, levels = c("AB", "BB", "AA"))
het <- ggplot() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 8, colour = "black"),
          axis.text.x = element_blank(),
          legend.position = "none",
          strip.text = element_blank(),
          strip.background = element_blank()) +
    scale_x_continuous(breaks = seq(0, 100000000, 10000000), labels = seq(0, 100, 10)) +
    geom_area(data = inversion_region, aes(x = x, y = ymax), fill = "grey") +
    #geom_area(data = inversion_region_29, aes(x = x, y = ymax), fill = "grey") +
    geom_line(data = td_het, aes(x = pos, y = Ho, colour = K2), size = 0.3) +
    scale_color_manual(values = c("AB" = col[1], "BB" = col[2], "AA" = col[3])) +
    ylab("Obs. Het.") +
    #labs(title = "B") +
    #geom_vline(xintercept = len_chr8) +
    facet_grid(~ chr, scales = "free_x", shrink = F, space = "free_x")
het


# MAF
frq.aa <- read.table("./01_data/chr8_29_AA.frq", header = T)
frq.aa <- frq.aa[c(n8:1, (n8+1):(n8+n29)),]
frq.ab <- read.table("./01_data/chr8_29_AB.frq", header = T)
frq.ab <- frq.ab[c(n8:1, (n8+1):(n8+n29)),]
frq.bb <- read.table("./01_data/chr8_29_BB.frq", header = T)
frq.bb <- frq.bb[c(n8:1, (n8+1):(n8+n29)),]

df_maf <- data.frame(chr = map[,1],
                     pos = map[,4],
                     "AA" = frq.aa$MAF,
                     "AB" = frq.ab$MAF,
                     "BB" = frq.bb$MAF)
td_maf <- gather(df_maf, "Karyotype", "MAF", -chr, -pos, factor_key = T)
td_maf$K2 <- factor(td_maf$Karyotype, levels = c("AB", "BB", "AA"))

maf <- ggplot() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 8, colour = "black"), 
          axis.text.x = element_blank(),
          #strip.text = element_text(size = 10, face = "bold", colour = "white"),
          #strip.background = element_rect(colour = "black", fill = "black"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "none") +
    scale_x_continuous(breaks = seq(0, 100000000, 10000000), labels = seq(0, 100, 10)) +
    scale_y_continuous(breaks = seq(0, 0.5, 0.25)) +
    geom_area(data = inversion_region, aes(x = x, y = ymin), fill = "grey") +
    geom_line(data = td_maf, aes(x = pos, y = MAF, colour = K2), size = 0.3) +
    scale_color_manual(values = c("AB" = col[1], "BB" = col[2], "AA" = col[3])) +
    ylab("MAF") +
    labs(title = "C") +
    #geom_vline(xintercept = len_chr8)
    facet_grid(~ chr, scales = "free_x", shrink = F, space = "free_x")
maf


# FST
fst.aa.bb <- read.table("./01_data/CAST_2016_50K_inversion_AA_BB.fst", header = T)
fst.aa.bb <- fst.aa.bb[c(n8:1, (n8+1):(n8+n29)),]
fst.ab.bb <- read.table("./01_data/CAST_2016_50K_inversion_AB_BB.fst", header = T)
fst.ab.bb <- fst.ab.bb[c(n8:1, (n8+1):(n8+n29)),]
fst.aa.ab <- read.table("./01_data/CAST_2016_50K_inversion_AA_AB.fst", header = T)
fst.aa.ab <- fst.aa.ab[c(n8:1, (n8+1):(n8+n29)),]

df_fst <- data.frame(chr = map[,1],
                     pos = map[,4],
                     "AAvBB" = fst.aa.bb$FST,
                     "AAvAB" = fst.aa.ab$FST,
                     "ABvBB" = fst.ab.bb$FST)
td_fst <- gather(df_fst, "Karyotype", "FST", -chr, -pos, factor_key = T)

fst <- ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 8, colour = "black"), 
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()) +
    scale_x_continuous(breaks = seq(0, 100000000, 10000000), labels = seq(0, 100, 10)) +
    geom_area(data = inversion_region, aes(x = x, y = ymax), fill = "grey") +
    geom_line(data = td_fst, aes(x = pos, y = FST, colour = Karyotype), size = 0.3) +
    scale_color_manual(values = c("AAvBB" = col[1], "AAvAB" = col[2], "ABvBB" = col[3])) +
    ylab("FST") +
    xlab("Fused chromosome position (Mb)") +
    #labs(title = "D") +
    #geom_vline(xintercept = len_chr8)
    facet_grid(~ chr, scales = "free_x", shrink = F, space = "free_x")
fst

#### put it all together ####
library(gridExtra)
pdf("./02_results/99_fusion/Figure2_fusion_wo_LD.pdf", width = 7, height = 2.5)
grid.arrange(grobs = list(pca, maf, het, fst), widths = c(1,2), heights = c(1.5,1,1.5), layout_matrix = rbind(c(1,2), c(1,3), c(1,4)))
dev.off()

