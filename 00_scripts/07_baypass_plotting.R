#### baypass ####

setwd("~/Desktop/CAST_2016_50K_analysis/")

library(ape)
library(corrplot)
source("/Applications/baypass_2.1/utils/baypass_utils.R")

omega <- as.matrix(read.table("./02_results/07_outliers/baypass/baypass_out_mat_omega.out"))
pop.names <- c("BAR","BHB","CAN","CWB","DUN","LNW","LSW","MSW","NWM","RBR","REN","SEV","TAX","UNW","USW")
dimnames(omega) <- list(pop.names, pop.names)
sample.order <- c(14,12,6,7,9,11,5,3,1,8,13,10,4,2,15)
omega <- omega[sample.order, sample.order]

# Visualize correlation matrix
cor.mat <- cov2cor(omega)
corrplot(cor.mat, method="color", mar=c(2,1,2,2)+0.1, main = expression("Correlation map based on"~hat(Omega)))

# Visualize correlation matrix as a tree
tree <- as.phylo(hclust(as.dist(1 - cor.mat**2)))
plot(tree, type = "p", main = expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

# Estimates of the XtX differentiation measure
baypass.snp.res <- read.table("./02_results/07_outliers/baypass/baypass_out_summary_pi_xtx.out", header = T)
plot(baypass.snp.res$M_XtX, pch = 16, cex = 0.3)

# Write the order of markers based on XtX
map <- read.table("./01_data/CAST_2016_50K.bim", stringsAsFactors =F)
write.table(map[baypass.snp.res$MRK[order(baypass.snp.res$M_XtX, decreasing = T)], 2], "./02_results/07_outliers/baypass_XtX_order.accnos", row.names= F, col.names = F, quote = F)




#### Simulate data to assess significance of outliers ####
pi.beta.coef <- read.table("./02_results/07_outliers/baypass/baypass_out_summary_beta_params.out", header = T)$Mean

# upload the original data to obtain total allele count
baypass.data <- geno2YN("./01_data/CAST_2016_50K.baypass")

# create the POD
setwd("./02_results/07_outliers/baypass/")
simu.baypass <- simulate.baypass(omega.mat = omega, nsnp = 10000, sample.size = baypass.data$NN, beta.pi = pi.beta.coef, pi.maf = 0, suffix = "baypass_pods")
setwd("../../../")
system("/Applications/baypass_2.1/sources/g_baypass -npop 15 -gfile ./02_results/07_outliers/baypass/G.baypass_pods -outprefix ./02_results/07_outliers/baypass/baypass_pods")

# read POD data and sanity check
pod.omega <- as.matrix(read.table("./02_results/07_outliers/baypass/baypass_pods_mat_omega.out"))
plot(pod.omega, omega)
abline(a = 0, b = 1)
fmd.dist(pod.omega, omega)

pod.pi.beta.coef <- read.table("./02_results/07_outliers/baypass/baypass_pods_summary_beta_params.out", header = T)$Mean
plot(pod.pi.beta.coef, pi.beta.coef)
abline(a = 0, b = 1)

# read POD XtX
pod.xtx <- read.table("./02_results/07_outliers/baypass/baypass_pods_summary_pi_xtx.out", header = T)$M_XtX
pod.thresh <- quantile(pod.xtx, probs = 0.99)

# Plot XtX with 1% threshold
plot(baypass.snp.res$M_XtX, pch = 16, cex = 0.3, col = ifelse((map[,1] %% 2) == 1, "black", "grey50"))
abline(h = pod.thresh, lwd = 2, lty = 2)

write.table(map[which(baypass.snp.res$M_XtX > pod.thresh), 2], "./02_results/07_outliers/baypass_outliers.accnos", row.names = F, col.names = F, quote = F)

