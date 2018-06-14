#### Spatial IBD-type analyses
library(vegan)
library(packfor)
library(reshape2)
library(adespatial)
library(robust)
library(MASS)
library(qvalue)
library(dplyr)
library(psych)
library(RColorBrewer)

setwd("~/Desktop/CAST_2016_50K_analysis/")

if (!dir.exists("./02_results/06_RDA/")) {
    dir.create("./02_results/06_RDA/")
}

# Read linear water distance among sites calculated using the riverdist package
geo_dist <- read.table("spatial_data/Ssa_mir_linear_river_distance_km.txt", header = T)
geo_dist <- as.matrix(geo_dist)

# Fix names of geographic distances
rownames(geo_dist) <- colnames(geo_dist) <- c("BAR", "BHB", "CAN", "CWB", "DUN", "LNW", "LSW", "MSW", "NWM", "RBR", "REN", "SEV", "TAX", "UNW", "USW")

# Create dbMEMs to explain spatial relationships
dbmem_coords <- dbmem(as.dist(geo_dist), MEM.autocor = "positive")

#### Enviro RDA ####
env_all <- read.csv("spatial_data/Enviro_data_extraction.csv")
row.names(env_all) <- as.character(env_all$DFO_siteID)

#### Prep data for sub-basin level RDA ####
env_sum <- env_all %>%
    group_by(Sub_Basin) %>%
    summarize(
        ELEV = weighted.mean(ELEV, Nind),
        DIST = weighted.mean(DIST, Nind),
        BC01 = weighted.mean(BC01, Nind),
        BC02 = weighted.mean(BC02, Nind),
        BC03 = weighted.mean(BC03, Nind),
        BC04 = weighted.mean(BC04, Nind),
        BC05 = weighted.mean(BC05, Nind),
        BC06 = weighted.mean(BC06, Nind),
        BC07 = weighted.mean(BC07, Nind),
        BC08 = weighted.mean(BC08, Nind),
        BC09 = weighted.mean(BC09, Nind),
        BC10 = weighted.mean(BC10, Nind),
        BC11 = weighted.mean(BC11, Nind),
        BC12 = weighted.mean(BC12, Nind),
        BC13 = weighted.mean(BC13, Nind),
        BC14 = weighted.mean(BC14, Nind),
        BC15 = weighted.mean(BC15, Nind),
        BC16 = weighted.mean(BC16, Nind),
        BC17 = weighted.mean(BC17, Nind),
        BC18 = weighted.mean(BC18, Nind),
        BC19 = weighted.mean(BC19, Nind)
    )
env_sum$LITH <- droplevels(env_all[env_all$rep == 1, "LITH", drop = T])
env <- as.data.frame(env_sum[,-1])
row.names(env) <- as.character(as.matrix(env_sum[,1]))

# read in all loci allele frequencies
all_frq <- read.table("./01_data/CAST_2016_50K.frq.strat", header = T)
all_freq <- dcast(all_frq, CLST ~ SNP, value.var = "MAF")
rownames(all_freq) <- as.character(all_freq[,1])
all_freq <- all_freq[rownames(env),-1]



# #### Prep data for sample site level analysis ####
# 
# # Need to generate allele frequencies for sampling sites
# # First create within file for calculating MAF with plink
# map <- read.table("./01_data/CAST_2016_50K.fam", header = F)
# map$Site <- sub("[A-Z].._", "", sub("_s.+", "", as.character(map[,2])))
# # List of sample sites to combine because of geographic proximity
# subs <- list(
#     LLSW_CatBr1 = "LLSW_CatBr2",
#     REN_no54 = "REN_rte108",
#     LSW_no61 = "LSW_Don",
#     RBR_rbr = "RBR_no92",
#     CAN_no74 = "CAN_Sabies",
#     MSW_no58 = "MSW_Porter",
#     TAX_no86 = "TAX_trib"
# )
# # Replace names for creating within file and calculate MAF
# for (l in names(subs)) {
#     map[map$Site == l, "Site"] <- subs[[l]]
# }
# write.table(map[,c(1,2,7)], "./01_data/CAST_2016_50K_site_within.txt", sep = " ", col.names = F, row.names = F, quote = F)
# system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --freq --within ./01_data/CAST_2016_50K_site_within.txt --out ./01_data/CAST_2016_50K_sample_site")
# 
# # Reduce env file to match
# env <- env_all[env_all$Nind > 12, ]
# env <- env[-which(as.character(env$DFO_siteID) %in% names(subs)), ]
# 
# 
# all_frq <- read.table("./01_data/CAST_2016_50K_sample_site.frq.strat", header = T)
# all_freq <- dcast(all_frq, CLST ~ SNP, value.var = "MAF")
# rownames(all_freq) <- as.character(all_freq[,1])
# all_freq <- all_freq[rownames(env),-1]


#### Common code for conducting enviro RDA ####

# # Create data type matrices
# temp <- as.matrix(env[, c("ELEV", "DIST", paste0("BC0", 1:9), "BC10", "BC11")])
# precip <- as.matrix(env[, paste0("BC", 12:19)])
# mig_dif <- scale(as.numeric(env[, "ELEV"] * env[, "DIST"] / 1000))
# lith <- env[, "LITH"]
# 
# # BioClim temperature variable PCA to reduce redundancy 
# temp_pca <- rda(temp, center = T, scale = T)
# temp_eig <- eigenvals(temp_pca)
# temp_keep <- sum(temp_eig > bstick(length(temp_eig), sum(temp_eig)))
# biplot(temp_pca, type = c("t","p"))
# 
# # BioClim precipitation variable PCA to reduce redundancy
# precip_pca <- rda(precip, center = T, scale = T)
# precip_eig <- eigenvals(precip_pca)
# precip_keep <- sum(precip_eig > bstick(length(precip_eig), sum(precip_eig)))
# biplot(precip_pca, type = c("t","p"))
# 
# # Create predictor dataframe
# X <- data.frame(temp = scores(temp_pca, choices = 1:temp_keep, display = "wa"),
#                 precip = scores(precip_pca, choices = 1:precip_keep, display = "wa"),
#                 #migd = mig_dif, # high VIF, correlated with MEM1 and temp PC1
#                 lith = lith)

X <- read.table("./01_data/summary_var_env_subbasin.txt")
row.names(X) <- X[,1]
X <- X[,-1]

pairs.panels(X)
loadings <- cor(env[,-22], X)
loadings[1:13, 2:3] <- NA
loadings[14:21, 1] <- NA
library(tidyverse)
library(ggplot2)
loadings_td <- gather(data.frame(env = rownames(loadings), loadings), PC, cor, -env, factor_key = T)
loadings_td$env <- factor(loadings_td$env, levels = rev(c("DIST","ELEV","BC01","BC02","BC03","BC04","BC05","BC06","BC07","BC08","BC09","BC10","BC11","BC12","BC13","BC14","BC15","BC16","BC17","BC18","BC19")))
loadings_td$PC <- factor(loadings_td$PC, levels = c("temp1_geo","prec1","prec2"))
loading_plot <- ggplot(loadings_td, aes(x = PC, y = env, fill = cor)) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(size = 8)) +
    geom_tile() +
    geom_text(aes(label = round(cor, 2))) +
    scale_fill_gradient2(low = "#f1a340", high = "#998ec3", limits = c(-1, 1), na.value = "grey") +
    scale_x_discrete("Environmental Summary Variables", labels = c("Temp./Elev.", "Winter precip.", "Summer precip."), expand = c(0,0)) +
    scale_y_discrete("Environmental Variables", expand = c(0,0))
ggsave("./02_results/Env_PCA_loadings.pdf", loading_plot, device = "pdf", width = 80, height = 120, units = "mm")
#env_rda <- rda(all_freq ~ PC1 + precip.PC1 + precip.PC2 + lith, X, scale = F) # ~1600 sec
env_rda <- rda(all_freq ~ ., X, scale = F) # ~1600 sec
vif.cca(env_rda)

# can parallelize this step, on a single processor (2.2GHz) takes ~26000 sec
env_stp_rda <- ordiR2step(rda(all_freq ~ 1, X, scale = F), scope = env_rda, Pin = 0.05, parallel = parallel::detectCores())
env_stp_rda2 <- ordistep(rda(all_freq ~ 1, X, scale = F), scope = env_rda, direction = "forward", Pin = 0.05, parallel = parallel::detectCores())

# This takes a long time so save output of env_stp_rda
saveRDS(env_stp_rda, "./02_results/07_outliers/env_stp_rda")

# Read env_stp_rda from RDS file
#env_stp_rda <- readRDS("./02_results/07_outliers/env_stp_rda")

#test overall significance and per axis
env_overall_test <- anova.cca(env_stp_rda, parallel = parallel::detectCores()) # 0.006
env_axis_test <- anova.cca(env_stp_rda, by="axis", parallel = parallel::detectCores()) # 0.001
env_margin_test <- anova.cca(env_stp_rda, by="margin", parallel = parallel::detectCores()) # 0.007
env_R2adj <- RsquareAdj(env_stp_rda) # r2 = 0.09471, r2.adj = 0.02507

# Print results
capture.output(summary(env_stp_rda), file = "./02_results/06_RDA/spatial_summary.txt")
capture.output(print(env_overall_test), file = "./02_results/06_RDA/spatial_overall_test.txt")
capture.output(print(env_axis_test), file = "./02_results/06_RDA/spatial_axis_test.txt")
capture.output(print(env_margin_test), file = "./02_results/06_RDA/spatial_margin_test.txt")
capture.output(print(env_R2adj), file = "./02_results/06_RDA/spatial_R2adj.txt")

# Create a plot
pal <- c(brewer.pal(12,"Paired"), "#000000", "#FFFFFF", "#CCCCCC")
pop <- factor(row.names(geo_dist))

pdf("./02_results/06_IBD/RDA_plot.pdf")
plot(env_stp_rda, scaling = 3, type = "n")
points(env_stp_rda, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)           # the SNPs
points(env_stp_rda, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = pal[pop]) # the sites
text(env_stp_rda, scaling = 3, display = "bp", col = "#0868ac", cex = 1)                           # the predictors
legend("bottomright", legend = levels(pop), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = pal)
dev.off()



#### Extract locus scores and calculate Mahalanobis dist etc following pcadapt ####

# extract raw scores
zscores <- scores(env_stp_rda, display = "species")

constr_axes <- grep("^RDA.",colnames(zscores))

# read SNP annotation to reorder zscores
snp_annot <- read.table("./01_data/CAST_2016_50K.map")
rownames(snp_annot) <- as.character(snp_annot$V2)

# if necessary, reorder SNPs to be in genomic order
if (!identical(as.character(snp_annot$V2), rownames(zscores))) {
    zscores <- zscores[as.character(snp_annot[,2]), ]
}

# center and scale scores
resscale <- apply(as.matrix(zscores[, constr_axes]), 2, scale)

qvals <- matrix(nrow = nrow(resscale), ncol = ncol(resscale), dimnames = list(rownames(zscores), colnames(zscores)[constr_axes]))

# Calculate significance based on robust Mahalanobis distance for each variable independently
for(i in 1:ncol(resscale)) {
    rcov <- cov.rob(resscale[, i]) # robust covariance
    resmaha <- mahalanobis(as.matrix(resscale[, i]), rcov$center, rcov$cov) # test statistic
    lambda <- median(resmaha) / qchisq(0.5, df = 1)
    pvals <- pchisq(resmaha / lambda, 1, lower.tail = FALSE)
    qvals[, i] <- qvalue(pvals)$qvalues
}

# Alternative is to use the 3 standard deviation cutoff as in Forester
outliers <- function(x,z){
    lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
    x < lims[1] | x > lims[2]               # locus names in these tails
}

candidate_outliers <- matrix(nrow = nrow(resscale), ncol = ncol(resscale), dimnames = list(rownames(zscores), colnames(zscores)[constr_axes]))
for(i in 1:length(constr_axes)) {
    candidate_outliers[, i] <- outliers(zscores[,i], 3)
}


# Write results file for reference
write.table(data.frame(raw_scores = zscores, qval = qvals, row.names = rownames(zscores)), "./02_results/06_RDA/ENV_RDA_SNP_scores.txt", quote = F)

# Manhattan plot of env_stp_rda score p-values
pdf("./02_results/06_RDA/env_rda_manhattan.pdf", width = 7, height = 3)
plot(1:length(pvals),-log10(pvals), col = c("grey50", "black")[(as.numeric(snp_annot[,1])%%2)+1],
     pch = 16, cex = 0.5, ylab = expression("-log"[10]*"(pvalue)"), xlab = "Chromosome", xaxt = "n")
axis(1, at = cumsum(table(snp_annot[,1])) - table(snp_annot[,1])/2, labels = 1:29, cex.axis = 0.7)
qsig<-which(qvals<0.05)
points(x = (1:length(pvals))[qsig], y = -log10(pvals)[qsig], col = 4, pch = 1, cex = 0.5)
dev.off()

# Create BED file of outliers for association testing
env_rda_outliers <- rownames(zscores)[which(qvals < 0.05)]

env_rda_outliers_bed <- data.frame(chr = paste0("chr",snp_annot[env_rda_outliers, 1]),
                                   startbp = as.numeric(snp_annot[env_rda_outliers, 4] - 1),
                                   endbp = as.numeric(snp_annot[env_rda_outliers, 4]),
                                   snp = env_rda_outliers)

write.table(env_rda_outliers_bed, "./01_data/env_rda_outliers.bed", sep = "\t", quote = F, col.names = F, row.names = F)

## Plotting highlighting SNPs
col_snps <- rep("#f1eef6", ncol(all_freq))
col_snps[which(colnames(all_freq) %in% env_rda_outliers)] <- "#e31a1c"
empty <- col_snps
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")

pdf("./02_results/06_RDA/RDA_plot_sig_snps.pdf")
plot(env_stp_rda, type="n", scaling=3, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
points(env_stp_rda, display="species", pch=21, cex=1, col="gray32", bg=col_snps, scaling=3)
points(env_stp_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(env_stp_rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("PC1"), bty="n", col="gray32", pch=21, cex=1, pt.bg = "#e31a1c")
dev.off()

#### Enviro association using MEM1 as conditioning variable ####

env_rda_C <- rda(all_freq ~ PC1 + precip.PC1 + precip.PC2 + lith + Condition(dbmem_coords[, "MEM1"]), X, scale = F)
vif.cca(env_rda_C)

#env_stp_rda_C <- ordiR2step(rda(all_freq ~ 1 + Condition(dbmem_coords[, "MEM1"]), X, scale = F), scope = env_rda_C, parallel = parallel::detectCores())
env_stp_rda_C <- ordistep(rda(all_freq ~ 1 + Condition(dbmem_coords[, "MEM1"]), X, scale = F), scope = env_rda_C, Pin = 0.05, direction = "forward", parallel = parallel::detectCores())

if (!dir.exists("./02_results/06_pRDA/")) {
    dir.create("./02_results/06_pRDA/")
}

# This takes a long time so save output of env_stp_rda
saveRDS(env_stp_rda_C, "./02_results/06_pRDA/env_stp_rda_C")

# Read env_stp_rda from RDS file
#env_stp_rda_C <- readRDS("./02_results/07_outliers/env_stp_rda_C")

#test overall significance and per axis
env_C_overall_test <- anova.cca(env_stp_rda_C, parallel = parallel::detectCores()) # 0.006
env_C_axis_test <- anova.cca(env_stp_rda_C, by="axis", parallel = parallel::detectCores()) # 0.001
env_C_margin_test <- anova.cca(env_stp_rda_C, by="margin", parallel = parallel::detectCores()) # 0.007
env_C_R2adj <- RsquareAdj(env_stp_rda_C) # r2 = 0.09471, r2.adj = 0.02507

# Print results
capture.output(summary(env_stp_rda_C), file = "./02_results/06_pRDA/spatial_summary.txt")
capture.output(print(env_C_overall_test), file = "./02_results/06_pRDA/spatial_overall_test.txt")
capture.output(print(env_C_axis_test), file = "./02_results/06_pRDA/spatial_axis_test.txt")
capture.output(print(env_C_margin_test), file = "./02_results/06_pRDA/spatial_margin_test.txt")
capture.output(print(env_C_R2adj), file = "./02_results/06_pRDA/spatial_R2adj.txt")

# Create a plot
pal <- c(brewer.pal(12,"Paired"), "#000000", "#FFFFFF", "#CCCCCC")
pop <- factor(row.names(geo_dist))

pdf("./02_results/06_pRDA/RDA_plot.pdf")
plot(env_stp_rda_C, scaling = 3, type = "n")
points(env_stp_rda_C, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)           # the SNPs
points(env_stp_rda_C, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = pal[pop]) # the sites
text(env_stp_rda_C, scaling = 3, display = "bp", col = "#0868ac", cex = 1)                           # the predictors
legend("bottomright", legend = levels(pop), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = pal)
dev.off()

# extract raw scores
zscores_C <- scores(env_stp_rda_C, display = "species")

constr_axes_C <- grep("^RDA.",colnames(zscores_C))


# if necessary, reorder SNPs to be in genomic order
if (!identical(as.character(snp_annot$V2), rownames(zscores_C))) {
    zscores_C <- zscores_C[as.character(snp_annot[,2]), ]
}

# center and scale scores
resscale_C <- apply(as.matrix(zscores_C[, constr_axes_C]), 2, scale)

qvals_C <- matrix(nrow = nrow(resscale_C), ncol = ncol(resscale_C), dimnames = list(rownames(zscores_C), colnames(zscores_C)[constr_axes_C]))

# Calculate significance based on robust Mahalanobis distance for each variable independently
for(i in 1:ncol(resscale_C)) {
    rcov_C <- cov.rob(resscale_C[, i]) # robust covariance
    resmaha_C <- mahalanobis(as.matrix(resscale_C[, i]), rcov_C$center, rcov_C$cov) # test statistic
    lambda_C <- median(resmaha_C) / qchisq(0.5, df = 1)
    pvals_C <- pchisq(resmaha_C / lambda_C, 1, lower.tail = FALSE)
    qvals_C[, i] <- qvalue(pvals_C)$qvalues
}


candidate_outliers_C <- matrix(nrow = nrow(resscale_C), ncol = ncol(resscale_C), dimnames = list(rownames(zscores_C), colnames(zscores_C)[constr_axes_C]))
for(i in 1:length(constr_axes_C)) {
    candidate_outliers_C[, i] <- outliers(zscores_C[,i], 3)
}


# Write results file for reference
write.table(data.frame(raw_scores = zscores_C, qval = qvals_C, row.names = rownames(zscores_C)), "./02_results/06_pRDA/ENV_pRDA_SNP_scores.txt", quote = F)

# Manhattan plot of env_stp_rda score p-values
pdf("./02_results/06_pRDA/env_rda_manhattan.pdf", width = 7, height = 3)
plot(1:length(pvals_C),-log10(pvals_C), col = c("grey50", "black")[(as.numeric(snp_annot[,1])%%2)+1],
     pch = 16, cex = 0.5, ylab = expression("-log"[10]*"(pvalue)"), xlab = "Chromosome", xaxt = "n")
axis(1, at = cumsum(table(snp_annot[,1])) - table(snp_annot[,1])/2, labels = 1:29, cex.axis = 0.7)
qsig_C<-which(qvals_C<0.05)
points(x = (1:length(pvals_C))[qsig_C], y = -log10(pvals_C)[qsig_C], col = 4, pch = 1, cex = 0.5)
dev.off()

# Create BED file of outliers for association testing
env_rda_outliers_C <- rownames(zscores_C)[which(qvals_C < 0.05)]

env_rda_outliers_C_bed <- data.frame(chr = paste0("chr",snp_annot[env_rda_outliers_C, 1]),
                                   startbp = as.numeric(snp_annot[env_rda_outliers_C, 4] - 1),
                                   endbp = as.numeric(snp_annot[env_rda_outliers_C, 4]),
                                   snp = env_rda_outliers_C)

write.table(env_rda_outliers_C_bed, "./01_data/env_rda_C_outliers.bed", sep = "\t", quote = F, col.names = F, row.names = F)


## Plotting highlighting SNPs
col_snps <- rep("#f1eef6", ncol(all_freq))
col_snps[which(colnames(all_freq) %in% env_rda_outliers_C)] <- "#e31a1c"
empty <- col_snps
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")

pdf("./02_results/06_pRDA/RDA_plot_sig_snps.pdf")
plot(env_stp_rda_C, type="n", scaling=3, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
points(env_stp_rda_C, display="species", pch=21, cex=1, col="gray32", bg=col_snps, scaling=3)
points(env_stp_rda_C, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(env_stp_rda_C, scaling=3, display="bp", col="#0868ac", cex=1)
#text(env_stp_rda_C, scaling=3, display="sites", col= brewer.pal(9, "Reds")[-1:-4][cut(X$PC1,5)], cex=1)
legend("bottomright", legend=c("PC1"), bty="n", col="gray32", pch=21, cex=1, pt.bg = "#e31a1c")
dev.off()