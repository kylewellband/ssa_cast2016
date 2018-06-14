#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
#args<-c("~/Desktop/CAST_2016_50K_analysis/02_results/CAST_2016_50K_neutral.res/CAST_2016_50K_neutral.xml", "CAST_2016_50K")
library(XML)
library(vegan)

# get PairFstMat adapted from R accessory file from Arlequin 3.5
getPWmat <- function(xmlfile) {
  tagData <- as.character(xmlValue(getNodeSet(xmlTreeParse(xmlfile, useInternal = T), "//PairFstMat")[[1]]))
  tagData2 <- strsplit(tagData, "\n")
  tagMat <- as.matrix(as.data.frame(tagData2))
  tagMat <- subset(tagMat, tagMat[,1] != "")
  tagMat <- tagMat[3:nrow(tagMat)]
  tagMat  <- gsub(" + ", " ", tagMat)
  Data <- strsplit(tagMat, " ")
  Row <- length(Data)
  Matrix <- as.matrix(as.data.frame(Data[1]))
  Matrix <- subset(Matrix, Matrix[,1] != "")
  Matrix <- rbind(Matrix, matrix(NA, ncol=1, nrow=(Row-1)))
  Matrix <- Matrix[2:nrow(Matrix),]
  numericList <- as.numeric(Matrix)
  numericMatrix <- t(as.matrix(numericList))
  if(Row >= 2){
    for(n in 2:(Row)){
      nextrow <- as.matrix(as.data.frame(Data[n]))
      nextrow <- subset(nextrow, nextrow[,1] != "")
      nextrow <- rbind(nextrow, matrix(NA, ncol=1, nrow=(Row-n)))
      nextrow <- nextrow[2:nrow(nextrow),]
      numericList <- as.numeric(nextrow)
      numericMatrix <- rbind(numericMatrix, t(as.matrix(numericList)))
    }
  }
  return(as.dist(numericMatrix))
}

getPWsig <- function(xmlfile) {
  tagData <- as.character(xmlValue(getNodeSet(xmlTreeParse(xmlfile, useInternal = T), "//PairFstPvalMat")[[1]]))
  tagData2 <- strsplit(tagData, "\n")
  tagMat <- as.matrix(as.data.frame(tagData2))
  tagMat <- subset(tagMat, tagMat[,1] != "")
  tagMat <- tagMat[2:nrow(tagMat)]
  tagMat  <- gsub("\\+-0.[0-9]+", " ", tagMat)
  tagMat  <- gsub("\\*", "0", tagMat)
  Data <- strsplit(tagMat, " ")
  Row <- length(Data)
  Matrix <- as.matrix(as.data.frame(Data[1]))
  Matrix <- subset(Matrix, Matrix[,1] != "")
  Matrix <- rbind(Matrix, matrix(NA, ncol=1, nrow=(Row-1)))
  Matrix <- Matrix[2:nrow(Matrix),]
  numericList <- as.numeric(Matrix)
  numericMatrix <- t(as.matrix(numericList))
  if(Row >= 2){
    for(n in 2:(Row)){
      nextrow <- as.matrix(as.data.frame(Data[n]))
      nextrow <- subset(nextrow, nextrow[,1] != "")
      nextrow <- rbind(nextrow, matrix(NA, ncol=1, nrow=(Row-n)))
      nextrow <- nextrow[2:nrow(nextrow),]
      numericList <- as.numeric(nextrow)
      numericMatrix <- rbind(numericMatrix, t(as.matrix(numericList)))
    }
  }
  return(as.dist(numericMatrix))
}

getPWnames <- function(xmlfile) {
  tagData <- as.character(xmlValue(getNodeSet(xmlTreeParse(xmlfile, useInternal = T), "//pairDistPopLabels")[[1]]))
  tagData2 <- strsplit(tagData, "\n")
  tagMat <- as.matrix(as.data.frame(tagData2))
  tagMat <- subset(tagMat, tagMat[,1] != "")
  tagMat <- tagMat[3:nrow(tagMat)]
  Data <- strsplit(tagMat, "\t")
  Matrix <- as.character(as.matrix(as.data.frame(Data)[2,]))
  return(Matrix)
}

# Read matrix from arlequin output
all_sample_pairwise_fst <- as.matrix(getPWmat(args[1]))
rownames(all_sample_pairwise_fst) <- colnames(all_sample_pairwise_fst) <- getPWnames(args[1])
npops <- nrow(all_sample_pairwise_fst)
all_sample_pairwise_fst_sig <- as.matrix(getPWsig(args[1]))
sample_order <- c(14,12,6,7,9,11,5,3,1,8,13,10,4,2,15)
all_sample_pairwise_fst <- all_sample_pairwise_fst[sample_order, sample_order]
write.table(all_sample_pairwise_fst, paste(args[2], ("_pwFst_heatmap.txt"), sep=""), quote = F, sep = "\t")
all_sample_pairwise_fst[all_sample_pairwise_fst_sig > 0.05/length(getPWsig(args[1]))] <- NA


pw_fst_heatmap <- function(Matrix = NULL, outfile = "") {
  
  if(is.null(Matrix) | !is.matrix(Matrix)) {stop("You must provide a valid FST matrix")}
  
  if(is.null(row.names(Matrix))) {
    Labels = NULL
  } else {
    Labels = row.names(Matrix)
  }
  
  mirror.matrix <- function(x) {
    xx <- as.data.frame(x);
    xx <- rev(xx);
    xx <- as.matrix(xx);
    xx;
  }
  
  #----Rotate matrix 270 clockworks----
  rotate270.matrix <- function(x) {
    mirror.matrix(t(x))
  }
  
  Matrix[upper.tri(Matrix)] = 0
  Matrix <- rotate270.matrix(Matrix)
  
  ColorRamp <- colorRampPalette(c("white", "steelblue1", "blue3"))
  
  #outfileGraphic <- paste(outfile, "pairFstMatrix ", timeAttr, ".png", sep="")
  outfileGraphic <- paste(outfile, ("_pwFst_heatmap.pdf"), sep="")
  
  #save graphic
  #png(outfileGraphic, width=1300, height=1300, res=144)
  pdf(outfileGraphic, width = 8, height = 8)  
  
  smallplot <- c(0.874, 0.9, 0.18, 0.83)
  bigplot <- c(0.13, 0.85, 0.14, 0.87)
  
  old.par <- par(no.readonly = TRUE)
  
  # draw legend --------------------------------
  par(plt = smallplot)
  
  # get legend values
  Min <- min(Matrix, na.rm=TRUE)
  Max <- max(Matrix, na.rm=TRUE)
  binwidth <- (Max - Min) / 64
  y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
  z <- matrix(y, nrow = 1, ncol = length(y))
  
  image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)
  
  # adjust axis if only one value exists
  if(Min == Max){
    axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
  } else {
    axis(side=4, las = 2, cex.axis=0.8)
  }
  
  box()
  mtext(text=expression(bold(F[ST])), side=4, line=2.5, cex=1.1)
  
  
  #draw main graphic ---------------------------
  a <- ncol(Matrix)
  b <- nrow(Matrix)
  
  x <- c(1:a)
  y <- c(1:b)
  
  par(new = TRUE, plt = bigplot)
  
  image(x,y,Matrix, col=ColorRamp(64),
        main=expression(bold(Matrix~of~pairwise~F[ST])), xlab="",
        ylab="", axes=FALSE)
  mmat <- ifelse(is.na(Matrix), 1, NA)
  image(x, y, mmat, axes = FALSE, xlab = "", ylab = "", 
        col = "grey85", add = TRUE)
  
  box()
  
  #add labels
  if(is.null(Labels)){
    axis(1, at = c(1:a))
    axis(2, at = c(1:b), labels=c(b:1))
    mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
          font=2)
    mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
          font=2)
  } else{
    axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
         las=2)
    axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
         las=2)
  }
  par(old.par)  #reset graphic parameters
  dev.off()
  
}

pw_fst_heatmap(as.matrix(all_sample_pairwise_fst), outfile = args[2])
