#### inveRsion Analysis ####
# This script uses inveRsion to identify the inversion break points in the fused 8 and 29 chromosomes
# it then uses these break points to classify the inversion status of individuals using invClust

# by: Kyle Wellband

setwd("~/Desktop/CAST_2016_50K_analysis/")

system("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --chr 8,29 --r2 square --recode A --make-bed --out ./01_data/chr8_29")

library(inveRsion)
library(snpStats)

data <- read.table("./01_data/chr8_29.raw", header = T, stringsAsFactors = F)
map <- read.table("./01_data/chr8_29.bim", header = F, stringsAsFactors = F)
n8 <- sum(map[, 1] == 8)
n29 <- sum(map[, 1] == 29)

geno <- as.matrix(data[, c(n8:1, (n8+1):(n8+n29))+6])+1

geno[is.na(geno)] <- 0

anno <- map[c(n8:1,(n8+1):(n8+n29)), ]
len_chr8 <- 26425719
anno[1:n8, 4] <- len_chr8 - anno[1:n8, 4]
anno[(n8+1):nrow(anno), 4] <- anno[(n8+1):nrow(anno), 4] + len_chr8
anno[,1] <- 8

# inveRsion

genodata_fusion <- setUpGenoDatSNPmat(Chr = 8, Geno = new("SnpMatrix", geno), Annot = anno, saveRes = F, saveGeno = T)

haploCode <- codeHaplo(genodata_fusion, blockSize = 5, minAllele = 0.1, saveRes = F)

scanRes <- scanInv(haploCode , window = 1, geno = T)

plot(scanRes)

invList <- listInv(scanRes, haploCode, geno = T, thBic = 0, saveRes = F, saveROI = F, all = T)

plot(invList, wROI = 1)

# we don't use this... its just for comparing the classification below
inv.genotype <- getClassif(invList, thBic = 500, wROI = 1, bin = T, geno = T, id = data$IID)


# invClust

library(invClust)

# define ROIs
# fusion is: 19175644 - 29484452 on inverted chr8 fused to chr29, or ends at 7250075 on normal chr8 and ends at 3058733 on normal chr29
# chr29 inversion is: 30235583 - 34092543 on fused chr, or 3809864 - 7666824 on chr29
roi<-data.frame(chr = c(8, 8), 
                LBP = c(invList@results[[1]]@intLeftBP[1] * 1E+6, invList@results[[2]]@intLeftBP[1] * 1E+6), 
                RBP = c(invList@results[[1]]@intRightBP[2] * 1E+6, invList@results[[2]]@intRightBP[2] * 1E+6), 
                reg = c("inv1", "inv2"))

names(anno) <- c("chromosome", "names", "cM", "position", "a1", "a2")

invcall.fus <- invClust(roi = roi, wh = 1, geno = new("SnpMatrix", geno), annot = anno, dim = 2)
invcall.29inv <- invClust(roi = roi, wh = 2, geno = new("SnpMatrix", geno), annot = anno, dim = 2)

plot(invcall.fus)
plot(invcall.29inv)

# Inversion status call
invclust.gt <- invGenotypes(invcall.fus)
levels(invclust.gt) <- c("AA", "AB", "BB")

invclust.gt.chr29 <- invGenotypes(invcall.29inv)
levels(invclust.gt.chr29) <- c("AA", "AB", "BB")

# comparing outputs
inv.gt.vec <- t(t(matrix(0:2)) %*% t(as.matrix(inv.genotype[,2:4])))
table(inv.gt.vec, invclust.gt)
accBic(invList, mem = as.numeric(invclust.gt)-1, nsub = 2*length(invclust.gt), npoints = 10, geno = T, wROI = 1)

# Write out some files for further analyses
write.table(cbind(data[,c("FID", "IID")], invclust.gt), "./01_data/CAST_2016_50K_inversion_status.txt", quote = F, col.names = F, row.names = F)

write.table(cbind(data[,c("FID", "IID")], invclust.gt.chr29), "./01_data/CAST_2016_50K_chr29_inversion_status.txt", quote = F, col.names = F, row.names = F)

write.table(table(data$FID, invclust.gt), "./01_data/Fusion_inversion_by_pop.txt", quote = F, col.names = T, row.names = T)

