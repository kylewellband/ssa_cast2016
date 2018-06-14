#### Fusion vs. PC2

library(ggplot2)

setwd("~/Desktop/CAST_2016_50K_analysis/")

fus_gt <- read.table("./01_data/CAST_2016_50K_inversion_status.txt")
fus_gt$fA <- numeric(length = dim(fus_gt)[1])
fus_gt$fA[fus_gt$V3 == "AA"] <- 2
fus_gt$fA[fus_gt$V3 == "AB"] <- 1

tab <- table(fus_gt$fA, fus_gt$V1)

trib_freq <- (tab[3,] * 2 + tab[2,]) / (2*colSums(tab))

env <- read.table("./01_data/summary_var_env_subbasin.txt")
names(env) <- c("Pop", "Temp./Elev.", "PC1", "PC2")

identical(as.character(env$Pop), names(trib_freq))

fus <- data.frame(env, fA = trib_freq, branch = factor(c(1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 1), labels = c("SW", "NW")))

pdf("./02_results/99_fusion/fusion_freq_summer_precip.pdf", width = 3.5, height = 3)
ggplot(fus, aes(x = PC2, y = fA, colour = branch)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 12, colour = "black"), legend.key = element_blank()) +
    geom_smooth(method = "lm", se = F, linetype = 2, show.legend = F) +
    geom_point(size = 3, show.legend = F) +
    geom_text(aes(label = Pop), colour = "black", nudge_x = 0.15, nudge_y = 0.01, size = 2) +
    ylab("Fusion Frequency (fA)") +
    xlab("Summer Precipitation") +
    scale_color_manual(values = c("SW" = "grey", "NW" = "black")) +
    geom_smooth(method = "lm", se = F, colour = "black")
dev.off()

