#### Which genes in the fusion LD block? ####

setwd("~/Desktop/CAST_2016_50K_analysis")

library(topGO)

system("bedtools intersect -a ./Ssalar_dbSNP/Ssal_chr_recode.gff3 -b ./01_data/fusion.bed | \
       grep -e 'gene' | \
       cut -f1,9 | \
       sed -e 's/ID=.*Name=//g' -e 's/\\.t[0-9]*//g' -e $'s/;Alias=/\t/g' -e $'s/;Note=/\t/g' | \
       sort | uniq | tee ./01_data/fusion.info | \
       cut -f2  | sort | uniq > ./01_data/fusion.accnos")

go.anno <- read.table("Ssalar_dbSNP/ssa_go_assoc.txt", header = T, stringsAsFactors = F)

inv.accnos <- as.character(read.table("./01_data/fusion.accnos", stringsAsFactors = F)[,1])

go.anno[,1] <- gsub("\\.t[0-9]*", "", go.anno[,1])

outliers.w.go <- inv.accnos[inv.accnos %in% go.anno[,1]]

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
gene.list <- factor(as.integer(names(gene2GO.map) %in% outliers.w.go))
names(gene.list) <- names(gene2GO.map)

fusion.GOdata <- new("topGOdata",
                      description = "Genes in fusion of chr 8 & 29", ontology = "BP",
                      allGenes = gene.list,
                      nodeSize = 5,
                      annot = annFUN.gene2GO, gene2GO = gene2GO.map)

fusion.result.Fisher <- runTest(fusion.GOdata, algorithm = "weight01", statistic = "fisher")
geneData(fusion.result.Fisher)
hist(score(fusion.result.Fisher), 50, xlab = "p-values")

# fusion.classic.Fisher <- runTest(fusion.GOdata, algorithm = "classic", statistic = "fisher")
# geneData(fusion.classic.Fisher)
# hist(score(fusion.classic.Fisher), 50, xlab = "p-values")

write.table(GenTable(fusion.GOdata, fusion.result.Fisher, topNodes = length(fusion.result.Fisher@score)), 
            "./02_results/99_fusion/fusion_topGO.txt", quote = F, row.names = F, sep = "\t")

printGraph(fusion.GOdata, fusion.result.Fisher, firstSigNodes = 5, fn.prefix = "./02_results/99_fusion/fusion_top5GO", pdfSW = T, useInfo = "np")
