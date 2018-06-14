#### Genotype assignment ####

setwd("~/Desktop/CAST_2016_50K_analysis/")
inds <- read.table("./01_data/CAST_2016_50K.fam", header = F, stringsAsFactors = F)[,1:2]


#### Sub-basin analysis ####

if(!dir.exists("./01_data/geneclass")) dir.create("./01_data/geneclass")
if(!dir.exists("./02_results/05_geneclass")) dir.create("./02_results/05_geneclass")

for (iter in 1:10) {
    
    if (file.exists(paste0("./01_data/geneclass/training", iter))) next
    
    training <- character()
    
    for (p in unique(inds[,1])) {
        training <- append(training, sample(as.character(inds[which(inds[,1] == p), 2]), size = ceiling(nrow(inds[which(inds[,1] == p),])/2), replace = F))
    }
    
    write.table(inds[which(inds[,2] %in% training),], paste0("./01_data/geneclass/training", iter), sep = " ", col.names = F, row.names = F, quote = F)
    
    system(paste0("/Applications/plink_mac/plink --bfile ./01_data/CAST_2016_50K --chr-set 29 --keep ./01_data/geneclass/training", iter, " --fst --family --out ./01_data/geneclass/training", iter))

    fst <- read.table(paste0("./01_data/geneclass/training", iter, ".fst"), header = T)
    fst <- fst[order(fst$FST, decreasing = T), ]
    
    for (nsnps in c(500, 1000, 2000, 3000, 4000, 5000)) {
        write.table(as.character(fst$SNP)[1:nsnps], paste0("./01_data/geneclass/training", iter, "_", nsnps, ".snps"), row.names = F, col.names = F, quote = F)
        system(paste0("/Applications/plink_mac/plink --file ./01_data/CAST_2016_50K --chr-set 29 --recode --extract ./01_data/geneclass/training", iter, "_", nsnps,".snps --out ./01_data/geneclass/training", iter, "_", nsnps))
        system(paste0("java -jar /Applications/PGDSpider_2.1.1.2/PGDSpider2-cli.jar -inputfile ./01_data/geneclass/training", iter, "_", nsnps, ".ped -outputfile ./01_data/geneclass/training", iter, "_", nsnps, ".txt -spid ./05_Misc/ped2genepop.spid"))
        system(paste0("./00_scripts/genepop_header.sh \"training", iter, "_", nsnps, "\" ./01_data/geneclass/training", iter, "_", nsnps, ".txt"))
    }
    
}

## Go away and run geneclass and or genodive for genotype assignment...

system("for f in ./01_data/geneclass/*.gdv; do ./00_scripts/process_genodive_res.sh $f; done")

assignment_success <- data.frame()
for (iter in 1:10) {
    fbase <- paste0("./01_data/geneclass/training", iter)
    training <- read.table(fbase, stringsAsFactors = F)
    holdout <- setdiff(as.character(inds$V2), training)
    
    exclude_first_gen_mig <- F
    for (nsnps in c(500, 1000, 2000, 3000, 4000, 5000)) {
        exclude <- read.table(paste0("./01_data/geneclass/results_training", iter, "_", nsnps, ".excluded"), header = T)
        exclude$Individual <- gsub("_$", "", exclude$Individual)
        row.names(exclude) <- exclude$Individual
        exclude$Current <- gsub("_.*", "", exclude$Current)
        exclude$Inferred <- gsub("_.*", "", exclude$Inferred)
        if(exclude_first_gen_mig) {
            true_holdout <- setdiff(holdout, exclude$Individual[exclude$Current != exclude$Inferred])
        } else {
            true_holdout <- holdout
        }
        ext <- xtabs(~ Current, data = exclude[exclude$Current != exclude$Inferred,])
        res <- read.table(paste0("./01_data/geneclass/results_training", iter, "_", nsnps, ".results"), header = T)
        res$Individual <- gsub("_$", "", res$Individual)
        row.names(res) <- res$Individual
        res$Current <- gsub("_.*", "", res$Current)
        res$Inferred <- gsub("_.*", "", res$Inferred)
        xt <- xtabs(~ Current + Inferred, data = res[true_holdout,])
        tmp <- matrix(0, nrow = dim(xt)[1], ncol = dim(xt)[1])
        rownames(tmp) = colnames(tmp) = row.names(xt)
        tmp[,colnames(xt)] <- xt
        xt <- tmp
        etmp <- numeric(dim(xt)[1])
        names(etmp) <- row.names(xt)
        etmp[names(ext)] <- ext
        ext <- etmp
        assignment_success <- rbind(assignment_success,
                                    data.frame(Nsnps = factor(rep(nsnps, length(rownames(xt)))),
                                               iter = factor(rep(iter, length(rownames(xt)))),
                                               pop = names(diag(xt)),
                                               success = diag(xt),
                                               tot = rowSums(xt),
                                               fg_mig = ext,
                                               per_suc = diag(xt)/rowSums(xt)))
    }
    
}


boxplot(per_suc ~ Nsnps, data = assignment_success, ylab = "Proportion of holdout set successfully assigned", xlab = "Number of SNPs")
boxplot(per_suc ~ pop, data = assignment_success, ylab = "Proportion of holdout set successfully assigned", xlab = "Population", las = 3)

library(ggplot2)

pdf("./02_results/05_geneclass/assignment_results_bwplot.pdf", width = 8, height = 4)
ggplot(data = assignment_success, aes(y = per_suc, x = pop, fill = Nsnps)) +
    geom_boxplot(position = position_dodge(width = 0.9)) +
    scale_fill_grey() +
    theme_linedraw() +
    theme(panel.grid = element_blank()) +
    ylim(c(0,1)) +
    ylab("Assignment success") +
    xlab("Sub-basin population") +
    geom_abline(slope = 0, intercept = 0.5, lty = 2) +
    geom_abline(slope = 0, intercept = 0.8, lty = 2)
dev.off()


library(dplyr)
tidy <- tbl_df(assignment_success)
tidy %>% group_by(pop) %>% summarise(mean(per_suc), min(per_suc), max(per_suc))
tidy %>% group_by(pop) %>% summarise(mean(fg_mig), min(fg_mig), max(fg_mig))


## Bubble plot
ass_res <- xt
pop_order <- c(14,12,6,7,9,11,5,1,3,8,13,10,4,2,15)
ass_res <- ass_res[pop_order,pop_order]

plot(1:15, 1:15, type = "n", xaxt = "n", yaxt = "n", xlab = "Assigned population", ylab = "Sampled population")
abline(a = 16, b = -1)

for(i in 1:15) {
    points(x = 1:15, y = rep(i, 15), pch = 16, col = "grey", cex = (ass_res[(15:1)[i],] / sum(ass_res[(15:1)[i],]))*3)
}

pdf("02_results/05_geneclass/bubble_plot.pdf", width = 6, height = 6)
ade4::table.value(ass_res, row.labels = row.names(ass_res), col.labels = row.names(ass_res))
dev.off()
