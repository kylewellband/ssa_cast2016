# Code to conduct NJ tree of chord distance adn perhaps DAPC

library(adegenet)
library(hierfstat)
library(ape)

setwd("~/Desktop/CAST_2016_50K_analysis/")

data <- read.PLINK("01_data/CAST_2016_50K_neutral.raw")

# Dch and NJ tree plotting
x <- as.matrix(data)
x[x == 0] <- "0101"
x[x == 1] <- "0102"
x[x == 2] <- "0202"
x[is.na(x)] <- "0000"

z <- df2genind(x, ncode = 2, pop = pop(data), NA.char = "00", type = "codom")

dch <- genet.dist(z, method = "Dch")
dch <- as.matrix(dch)
rownames(dch) <- colnames(dch) <- popNames(z)

plot(nj(as.dist(dch)))
add.scale.bar()

# DAPC

# basic codes, no structure found... maybe K=2 but this is similar to Admixture...
grp <- find.clusters(data, max.n.clust = 10)
dapc <- dapc(x = data, pop = grp$grp)

sub <- character()
for(p in popNames(data)) {
  sub <- append(sub, sample(indNames(data)[which(pop(data) == p)], 5, replace = F))
}

data_sub <- data[(indNames(data) %in% sub),]
grp_sub <- pop(data_sub)
test <- dapc(x = data_sub, grp_sub)


# On all inds
dapc <- dapc(x = data, pop = pop(data)) # kept 200 PCs and 7 DFs
scatter(dapc, xax = 1, yax = 2)
scatter(dapc, xax = 1, yax = 3)
scatter(dapc, xax = 4, yax = 7)
scatter(dapc, xax = 1, yax = 5)
scatter(dapc, xax = 1, yax = 6)
scatter(dapc, xax = 1, yax = 7)

plot(density(dapc$var.contr[,1]))
loi <- intersect(which(dapc$var.contr[,1]>0.0004), which(dapc$var.contr[,2]>0.0004))
doi<-as.matrix(tapply(as.matrix(data[,loi[3]]), pop(data), function(i) mean(i, na.rm = T)/2))
doi <- cbind(doi, 1-doi)

matplot(doi, pch=c("a","c"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1.5, main="", ylim = c(0,1))
axis(1, 1:17, popNames(data), las = 2)

as.matrix(data[,intersect(which(dapc$var.contr[,1]>0.0004), which(dapc$var.contr[,2]>0.0004))])

sub_data<-data[,union(union(union(union(which(dapc$var.contr[,1]>0.0005), which(dapc$var.contr[,2]>0.0005)), which(dapc$var.contr[,3]>0.0005)), which(dapc$var.contr[,5]>0.0005)), which(dapc$var.contr[,6]>0.0005))]
grp<-find.clusters(data, max.n.clust = 10)

sub_dapc<-dapc(sub_data, grp$grp)

barplot(t(sub_dapc$posterior), col = 2:5, space = 0, border = NA)
