setwd("~/Desktop/CAST_2016_50K_analysis/")
fusion8_29 <- read.table("01_data/CAST_2016_50K_inversion_status.txt", stringsAsFactors = F)

fusion8_29[fusion8_29[,3] == "AA", 4] <- 2
fusion8_29[fusion8_29[,3] == "AB", 4] <- 1
fusion8_29[fusion8_29[,3] == "BB", 4] <- 0

CreatePlinkFileFromInversion <- function(inv.status = NULL, outfile.prefix = NULL, chr = NULL, pos = NULL, sex = NULL, pheno = NULL) {
    
    if (is.null(inv.status)) {
        stop("Error: You must provide an inversion status file")
    }
    if (is.null(chr)) {
        print("Warning: You did not provide a chromosome number. Chr has been assigned a value of 1")
        chr <- 1
    }
    if (is.null(pos)) {
        print("Warning: You did not provide a base pair postion number. It has been assigned a value of 1")
        pos <- 1
    }
    if (is.null(sex)) {
        sex <- numeric(nrow(inv.status))
    }
    if (is.null(pheno)) {
        pheno <- numeric(nrow(inv.status))
    }
    
    gt.vec <- inv.status[,3]
    
    if (is.factor(gt.vec)) {
        genotypes <- levels(gt.vec)
    } else {
        genotypes <- unique(gt.vec)
    }
    
    genotypes <- sort(genotypes)
    
    for (i in 1:length(genotypes)) {
        if (identical(unlist(strsplit(genotypes[i] ,NULL))[1], unlist(strsplit(genotypes[i] ,NULL))[2])) {
            next
        } else {
            het <- i
            break
        }
    }
    
    hom <- c(1:3)[!(1:3 %in% het)]
    
    ped.gt <- character(nrow(inv.status))
    ped.gt[gt.vec == genotypes[hom[1]]] <- "1 1"
    ped.gt[gt.vec == genotypes[het]] <- "1 2"
    ped.gt[gt.vec == genotypes[hom[2]]] <- "2 2"
    
    map <- cbind(chr, outfile.prefix, 0, pos)
    
    ped <- cbind(inv.status[,1:2], matrix(0, nrow = nrow(inv.status), ncol = 2), sex, pheno, ped.gt)
    
    if (is.null(outfile.prefix)) {
        print("Warning: You did not provide an outfile prefix. File will be saved in the working directory as: \"inversion.ped\"")
        outfile.prefix <- "inversion"
    }
    
    write.table(map, paste0(outfile.prefix, ".map"), sep = " ", col.names = F, row.names = F, quote = F)
    write.table(ped, paste0(outfile.prefix, ".ped"), sep = " ", col.names = F, row.names = F, quote = F)
    
}

CreatePlinkFileFromInversion(inv.status = fusion8_29, "./01_data/fusion8_29", chr = 8, pos = 1)
