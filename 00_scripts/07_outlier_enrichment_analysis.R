setwd("~/Desktop/CAST_2016_50K_analysis/")

library(topGO)

outliers <- read.table("./01_data/env_rda_C_outliers.accnos", stringsAsFactors = F)[,1]
all.50K.genes <- read.table("./Ssalar_dbSNP/CAST_50K_10000_window.accnos", stringsAsFactors = F)[,1]
go.anno <- read.table("ssa_go_assoc.txt", header = T, stringsAsFactors = F)
go.anno[,1] <- gsub("\\.t[0-9]*", "", go.anno[,1])

outliers.w.go <- outliers[outliers %in% go.anno[,1]]
all.50K.genes.w.go <- all.50K.genes[all.50K.genes %in% go.anno[,1]]

CreateGene2GOMapping <- function(x) {
    nr <- nrow(x)
    map <- list()
    for(r in 1:nr) {
        map[[as.character(x[r, 1])]] <- append(map[[as.character(x[r, 1])]], x[r,2])
    }
    return(map)
}

# create GO mapping, takes about 20 minutes
gene2GO.map <- CreateGene2GOMapping(go.anno)

gene2GO.map <- gene2GO.map[names(gene2GO.map) %in% all.50K.genes.w.go]
gene.list <- factor(as.integer(all.50K.genes.w.go %in% outliers.w.go))
names(gene.list) <- all.50K.genes.w.go

env.rda.GOdata <- new("topGOdata",
                    description = "env RDA outliers", ontology = "BP",
                    allGenes = gene.list,
                    nodeSize = 5,
                    annot = annFUN.gene2GO, gene2GO = gene2GO.map)

env.rda.result.Fisher <- runTest(env.rda.GOdata, algorithm = "weight01", statistic = "fisher")

geneData(env.rda.result.Fisher)
hist(score(env.rda.result.Fisher), 50, xlab = "p-values")

# env.rda.classic.Fisher <- runTest(env.rda.GOdata, algorithm = "classic", statistic = "fisher")
# geneData(env.rda.classic.Fisher)
# hist(score(env.rda.classic.Fisher), 50, xlab = "p-values")

write.table(GenTable(env.rda.GOdata, env.rda.result.Fisher, topNodes = length(env.rda.result.Fisher@score)), 
            "./02_results/07_outliers/env_rda_C_top5GO.txt", quote = F, row.names = F, sep = "\t")

printGraph(env.rda.GOdata, env.rda.result.Fisher, firstSigNodes = 5, fn.prefix = "./02_results/07_outliers/env_rda_C_top5GO", pdfSW = T, useInfo = "np")

