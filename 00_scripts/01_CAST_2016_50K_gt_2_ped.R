#!/usr/bin/env Rscript
# CAST_2016_50K_gt_2_ped.R
# by: Kyle Wellband
#
# This script was written to take as input the raw "*.txt" genotype output from 
# the Axiom Genotyping Console software by CIGENE, the SNP positions with offset
# from Cooke / Aquagen and sample annotation information from UNB / CAST partners
# to create PLINK format ("*.ped" and "*.map") data files for the analysis of
# populations genetic structure of Atlantic Salmon in the Miramichi River.
#

# Set working directory to location of all files
setwd("~/Desktop/CAST_2016_50K_analysis/01_Data/rawdata/")

# Read good calls at CR>0.97
x <- read.table("AxiomGT1.calls.txt", header = T, comment.char = "#")

# Read recovered calls at CR>0.9
y <- read.table("AxiomGT1.calls.mod_140bad_samples.txt", header = T, comment.char = "#")

# Check order of SNPs is the same
if (sum(as.character(x$probeset_id) != as.character(y$probeset_id))) {
  stop("The order of SNPs differs among sample files!")
}

# Combine two genotype files together
z <- cbind(x,y[,-1])
row.names(z) <- as.character(z$probeset_id)

# Read SNP position information and reorder to reflect genomic order
pos_info <- read.table("NA50k_vs_ICSASG_v2.positions.100bprandomoffset_w_adjustments.txt", header = F, row.names = 1, colClasses = c("character"))
pos_info <- pos_info[-which(pos_info[,1] == "-"),]
pos_info <- pos_info[order(as.numeric(pos_info[,2])),]
pos_info <- pos_info[order(pos_info[,1]),]
pos_info[,1] <- gsub(pattern = "ssa", replacement = "chr", x = pos_info[,1])

# Create "*.map" file with chromosomes and indexed genomic positions
map <- data.frame(chr = pos_info[,1], snp_id = row.names(pos_info), dist = rep(0, nrow(pos_info)), bp = pos_info[,2])
write.table(map, "../CAST_2016_50K_all.map", quote = F, row.names = F, col.names = F)

# Separate, reorder, and transpose genotype data
gt <- t(z[row.names(pos_info),-1])

# Read in file that links population (sub-basin) origin of fish as well as sexto "sample_XXXX" codes from genotyping
sample_anno <- read.table("sample_anno.txt", header = T, row.names = 1, stringsAsFactors = F)

# Convert letter sex IDs into binary format for PLINK
sample_anno$sex_bin <- as.numeric(sample_anno$sex == "F") + 1

# Reduce annotation to those samples successfully sequenced
sample_anno <- sample_anno[which(row.names(sample_anno) %in% row.names(gt)), ]

# Convert high quality (CR>0.97) and low quality (CR<0.9) genotypes into phenotype info
pheno <- rep(c(0,1), c(ncol(x)-1, ncol(y)-1)) + 1
names(pheno) <- row.names(gt)
pheno <- pheno[row.names(sample_anno)]

# Create new individual codes
ind <- paste(sample_anno$new_subbasin_assignment, 
             sample_anno$Site_code,
             gsub("_.*", "", gsub("Sample_", "s", row.names(sample_anno))),
             sep = "_")

# Reorder genotypes to match annotation
gt <- gt[row.names(sample_anno),]

# Recode single digit genotype codes into two allele genotype codes
gt[gt == 0] = "1 1" # reference homozygote
gt[gt == 1] = "1 2" # heterozygote
gt[gt == 2] = "2 2" # variant homozygote
gt[gt == -1] = "0 0" # missing

ped <- data.frame(pop = sample_anno$new_subbasin_assignment,
                  ind = ind,
                  sire = rep(0, dim(gt)[1]),
                  dam = rep(0, dim(gt)[1]),
                  sex = sample_anno[, "sex_bin"],
                  pheno = pheno,
                  gt)

ped <- ped[order(ped$pop, ped$ind), ]

write.table(ped, "../CAST_2016_50K_all.ped", quote = F, row.names = F, col.names = F)
