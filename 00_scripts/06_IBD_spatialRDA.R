#### Spatial IBD-type analyses
library(vegan)
library(reshape2)
library(adespatial)
library(robust)
library(MASS)
library(qvalue)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/CAST_2016_50K_analysis/")

# Read linear water distance among sites calculated using the riverdist package
geo.dist <- read.table("spatial_data/Ssa_mir_linear_river_distance_km.txt", header = T)
geo.dist <- as.matrix(geo.dist)

# Fix names of geographic distances
rownames(geo.dist) <- colnames(geo.dist) <- c("BAR", "BHB", "CAN", "CWB", "DUN", "LNW", "LSW", "MSW", "NWM", "RBR", "REN", "SEV", "TAX", "UNW", "USW")

# Read genetic distance (FST) and linearize for mantel test
genet.dist <- read.table("02_results/00_Basic_Stats/CAST_2016_50K_pwFst.txt")
genet.dist <- as.matrix(genet.dist)
genet.dist <- genet.dist[rownames(geo.dist), colnames(geo.dist)]
linear.genet.dist <- genet.dist/(1-genet.dist)

globalFST <- mean(genet.dist) # approximates global FST

# Perform mantel test
mantel_res <- mantel(as.dist(geo.dist), as.dist(linear.genet.dist))

capture.output(print(mantel_res), file = "./02_results/06_IBD/mantel_results.txt")

ibd_df <- data.frame(geo = as.numeric(as.dist(geo_dist)),
                     gen = as.numeric(as.dist(linear_genet_dist)))

pdf("./02_results/06_IBD/IBD_mantel.pdf", height = 3, width = 3.5)
ggplot(ibd_df, aes(x = geo, y = gen)) +
    geom_point() +
    ylab(expression("F"[ST]*" / (1 - F"[ST]*")")) +
    xlab("Linear river distance (KM)") +
    geom_smooth(color = "black", size = 1, method = "lm") +
    theme_classic() +
    theme(panel.grid = element_blank())
dev.off()


#### RDA of spatial effects on allele frequencies #### 

# spatial dbMEM procedure from Forester et al. 2017. bioRxiv

# If necessary use PLINK to calculate population allele frequencies for lowlink dataset
if (!file.exists("./01_data/CAST_2016_50K_neutral.frq.strat")) {
    system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K_neutral --chr-set 29 --freq --family --out ./01_data/CAST_2016_50K_neutral")
}

# Read and reshape population allele frequencies
ll.frq <- read.table("./01_data/CAST_2016_50K_neutral.frq.strat", header = T)
ll.freq <- dcast(ll.frq, CLST ~ SNP, value.var = "MAF")
rownames(ll.freq) <- as.character(ll.freq[,1])
ll.freq <- ll.freq[rownames(geo.dist),-1]

# calculate dbMEM from matrix of linear river distances
dbmem.coords <- dbmem(as.dist(geo.dist), MEM.autocor = "positive")

#perform initial RDA using all MEMs
spatial.rda <- rda(ll.freq ~ ., data.frame(dbmem.coords), scale = F) 
vif.cca(spatial.rda)

# perform forward step selection of MEMs retaining those significant (p < 0.05)
# can parallelize this step if necessary
spatial.stp.rda <- ordistep(rda(ll.freq ~ 1, data.frame(dbmem.coords)), scope = formula(spatial.rda), scale = F, direction = "forward", parallel = parallel::detectCores())
spatial.selected <- attributes(spatial.stp.rda$terms)$term.labels

#test overall significance and per axis
spatial.stp.rda <- rda(ll.freq ~ MEM1, data.frame(dbmem.coords), scale = F) 
spatial.overall.test <- anova.cca(spatial.stp.rda, parallel = parallel::detectCores()) # 0.031
spatial.axis.test <- anova.cca(spatial.stp.rda, by = "axis", parallel = parallel::detectCores()) # 0.033
spatial.margin.test <- anova.cca(spatial.stp.rda, by = "margin", parallel = parallel::detectCores()) # 0.032
spatial.R2adj <- RsquareAdj(spatial.stp.rda) # 0.015

capture.output(summary(spatial_stp_rda), file = "./02_results/06_IBD/spatial_summary.txt")
capture.output(print(spatial_overall_test), file = "./02_results/06_IBD/spatial_overall_test.txt")
capture.output(print(spatial_axis_test), file = "./02_results/06_IBD/spatial_axis_test.txt")
capture.output(print(spatial_margin_test), file = "./02_results/06_IBD/spatial_margin_test.txt")
capture.output(print(spatial_R2adj), file = "./02_results/06_IBD/spatial_R2adj.txt")


# Print RDA triplot
sc <- scores(spatial.stp.rda, scaling = 3, display = c("sp", "wa", "bp"))

pdf("./02_results/06_IBD/spatial_RDA_triplot.pdf", width = 3.5, height = 3.5)
ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 12, colour = "black"), legend.text = element_text(size = 12, colour = "black")) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_point(data = as.data.frame(sc$species), aes(x = RDA1, y = PC1), color = "grey32", size = 0.5, shape = "+") +
    geom_point(data = as.data.frame(sc$sites), mapping = aes(x = RDA1, y = PC1), size = 3, colour = c("grey", "black")[c(c(1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 1))]) +
    geom_text(data = as.data.frame(sc$sites), mapping = aes(x = RDA1, y = PC1, label = rownames(as.data.frame(sc$sites))), nudge_x = 0.1, nudge_y = 0.1, size = 2) +
    geom_segment(data = as.data.frame(sc$biplot), aes(x = 0, y = 0, xend = RDA1, yend = PC1), size = 1, arrow = arrow(length = unit(0.03, "npc"))) +
    geom_label(data = as.data.frame(sc$biplot), mapping = aes(x = RDA1, y = PC1, label = "MEM1"), size = 2, nudge_x = -0.1)
dev.off()

