setwd("~/Desktop/CAST_2016_50K_analysis/")

inv <- read.table("01_data/CAST_2016_50K_inversion_status.txt")
env <- read.csv("spatial_data/Enviro_data_extraction.csv")
env <- read.table("01_data/Fusion_by_site_env.txt", header = T)
rownames(env) <- as.character(env$DFO_siteID)
inv$DIST <- env[sub("_s.*", "", sub("[A-Z][A-Z][A-Z]_", "", as.character(inv[,2]))), "DIST"]

boxplot(inv$DIST ~ inv$V3)
plot(as.numeric(inv$V3), inv$DIST)

xtab<-table(sub("_s.*", "", sub("[A-Z][A-Z][A-Z]_", "", as.character(inv[,2]))), inv$V3)
xtab<-xtab[-which(rowSums(xtab)<10),]
freq<-((xtab[,1]*2) + xtab[,2]) / (2*rowSums(xtab))
indep.var <- env[names(freq),"BC18"] # BC15 & 16 & 18
lm1 <- lm(freq ~ indep.var, weights = rowSums(xtab))
summary(lm1)
plot(freq ~ indep.var, main = sprintf("R2 = %.3f, p = %.3f", summary(lm1)$adj.r.squared, anova(lm1)$'Pr(>F)'[1]))
abline(lm1)

# by sub-basin
env <- env[env$rep == 1,]
row.names(env) <- as.character(env$Sub_Basin)
xtab<-table(inv[,c(1,3)])
freq<-((xtab[,1]*2) + xtab[,2]) / (2*rowSums(xtab))

temp <- as.matrix(env[names(freq), c(paste0("BC0", 1:9), "BC10", "BC11")])
precip <- as.matrix(env[names(freq), paste0("BC", 12:19)])
elev <- as.numeric(env[names(freq), "ELEV"])
lith <- as.integer(env[names(freq), "GEO_KEY"])

# BioClim temperature variable PCA to reduce redundancy
temp.pca <- vegan::rda(temp, center = T, scale = T)
temp.eig <- vegan::eigenvals(temp.pca)
temp.keep <- sum(temp.eig > vegan::bstick(length(temp.eig), sum(temp.eig)))
biplot(temp.pca, type = c("t","p"))

# BioClim precipitation variable PCA to reduce redundancy
precip.pca <- vegan::rda(precip, center = T, scale = T)
precip.eig <- vegan::eigenvals(precip.pca)
precip.keep <- sum(precip.eig > vegan::bstick(length(precip.eig), sum(precip.eig)))
biplot(precip.pca, type = c("t","p"))

# Create predictor dataframe
X <- data.frame(temp = vegan::scores(temp.pca, choices = 1:temp.keep, display = "wa"),
                precip = vegan::scores(precip.pca, choices = 1:precip.keep, display = "wa"),
                #elev = elev,
                mig.dif = scale(as.numeric(env[names(freq), "ELEV"] * env[names(freq), "DIST"] / 1000)),
                lith = lith)

# for individual variables
indep.var <- env[names(freq),"BC18"] # BC15 & 16 & 18

# for PCs 
v = 5
indep.var <- X[,v]
lm1 <- lm(freq ~ indep.var, weights = rowSums(xtab))
plot(freq ~ indep.var, main = sprintf("R2 = %.3f, p = %.3f", summary(lm1)$adj.r.squared, anova(lm1)$'Pr(>F)'[1]), xlab = colnames(X)[v])
abline(lm1)


levels(inv$V1) <- c("SW", "SW", "SW", "SW", "SW", "NW", "NW", "SW", NA, "SW", "SW", "NW", "SW", "NW", "SW")
inv$fA <- 0
inv$fA[inv$V3 == "AA"] <- 2
inv$fA[inv$V3 == "AB"] <- 1

chisq.test(table(inv$fA, inv$V1))




#### Plotting for Allen
env_all <- read.csv("spatial_data/Enviro_data_extraction.csv")
row.names(env_all) <- as.character(env_all$DFO_siteID)

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

library(tidyverse)
library(ggplot2)
df <- data.frame(Sub_Basin = factor(rownames(env)), freq = freq, env[,-22])
td <- gather(df, key = "env", value = "val", -freq, -Sub_Basin, factor_key = T)

pdf("~/Desktop/fusion_by_each_env_variable.pdf", width = 7.5, height = 10)
pp <- ggplot(td, aes(x = val, y = freq)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black")) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    facet_wrap(~ env, ncol = 4, scales = "free_x") +
    ylab("Fusion frequency (fA)") +
    xlab("")
pp
dev.off()



env<-read.table("01_data/summary_var_env_subbasin.txt", header = T)
rownames(env) <- env[,1]
xtab<-table(sub("_.*", "", sub("[A-Z][A-Z][A-Z]_", "", as.character(inv[,2]))), inv$V3)
xtab<-xtab[-which(rowSums(xtab)<12),]
freq<-((xtab[,1]*2) + xtab[,2]) / (2*rowSums(xtab))

summary(lm(freq ~ env$prec2))
plot(env$prec2, freq, ylim = c(0,1), xlim = c(-2,2))


## Kyle's hacked logistic regression
linker <- sub("_s.*", "", sub("[A-Z][A-Z][A-Z]_", "", as.character(inv[,2])))

response <- numeric(nrow(inv))
response[inv[,3] == "AA" | inv[,3] == "AB"] <- 1
pred <- env[sub("_.*", "", as.character(inv[,2])), "prec2"]
m1 <- glm(response ~ pred, family = binomial(link = "logit"))
anova(m1, test = "Chisq")
library(pscl)
pR2(m1)
plot(pred, response)
