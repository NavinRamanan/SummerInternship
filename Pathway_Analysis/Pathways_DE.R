#Navin Ramanan pathway analysis code


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")
install.packages("BiocManager")



library(rWikiPathways)


installApp('WikiPathways') 
installApp('CyTargetLinker') 
installApp('stringApp') 
installApp('enrichmentMap')

load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggplot",
  "dotplot",
  "goplot",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}


library(rWikiPathways)
library(statmod)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(purrr)
library(tidyverse)
library(biomaRt)
library(dotplot)
library(goplot)

#pathway analysis
Cortex.expr <- read.csv("/Users/navinramanan/Desktop/SI/SI_Files/APP_all/Hippocampus_deseq_results_APP_all_padj_05_filtered.csv", stringsAsFactors = FALSE)

nrow(Cortex.expr)
head(Cortex.expr)

de.genes <- Cortex.expr[Cortex.expr$log2FoldChange > 1 | Cortex.expr$log2FoldChange < -1 & Cortex.expr$padj < 0.1, 1]

bkgd.genes <- Cortex.expr[,1]
view(de.genes)


de.genes.entrez <- clusterProfiler::bitr(de.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

bkgd.genes.entrez <- bitr(bkgd.genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

#gene ontology analysis for all DE genes
egobp <- clusterProfiler::enrichGO(
  gene     = de.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Mm.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)



head(egobp,10)

print(de.genes)

barplot(egobp, showCategory = 20)
dotplot(egobp, showCategory = 20)
goplot(egobp)
view(egobp)

## extract a dataframe with results from object of type enrichResult
egobp.results.df <- egobp@result
view(egobp.results.df)
## create a new column for term size from BgRatio
egobp.results.df$term.size <- gsub("/(\\d+)", "", egobp.results.df$BgRatio)
view(egobp.results.df$term.size)

## filter for term size to keep only term.size => 3, gene count >= 5 and subset
egobp.results.df <- egobp.results.df[which(egobp.results.df[,'term.size'] >= 3 & egobp.results.df[,'Count'] >= 5),]
egobp.results.df <- egobp.results.df[c("ID", "Description", "pvalue", "qvalue", "geneID")]

## format gene list column
egobp.results.df$geneID <- gsub("/", ",", egobp.results.df$geneID)

## add column for phenotype
egobp.results.df <- cbind(egobp.results.df, phenotype=1)
egobp.results.df <- egobp.results.df[, c(1, 2, 3, 4, 6, 5)]

## change column headers
colnames(egobp.results.df) <- c("Name","Description", "pvalue","qvalue","phenotype", "genes")

egobp.results.filename <-file.path(getwd(),paste("clusterprofiler_cluster_enr_results.txt",sep="_"))
write.table(egobp.results.df,egobp.results.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.1", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',egobp.results.filename ,
                   sep=" ")

#enrichment map command will return the suid of newly created network.

#GO analysis for all DE genes
ewp.de <- clusterProfiler::enrichWP(
  de.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  organism = "Mus musculus",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, 
)

ewp.de <- DOSE::setReadable(ewp.de, org.Mm.eg.db, keyType = "ENTREZID")
head(ewp.de)

barplot(ewp.de, showCategory = 20)
dotplot(ewp.de, showCategory = 20)

#Find pathways with specified name. Will output information in console and save pathway ids 
#findPathwayNamesByText("TyroBP")
#lc.pathways <- findPathwaysByText("TyroBP")  #quotes inside query to require both terms
#mmu.lc.pathways <- lc.pathways %>% 
 # dplyr::filter(species == "Mus musculus") # 
#mmu.lc.pathways$name # display the pathway titles
#lc.wpids <- mmu.lc.pathways$id
#lc.wpids

ewp.de.wpids <- ewp.de$ID
ewp.de.wpids
#Create pathway map for specified pathway. WP3625 is the TyroBP causal network in microglia pathway

#url <- getPathwayInfo("WP3625")$url
#browseURL(url)
#RCy3::commandsRun('wikipathways import-as-pathway id=WP3625') 
#toggleGraphicsDetails()
#loadTableData(Cortex.expr, data.key.column = "X", table.key.column = "Ensembl")
#setNodeColorMapping("log2FoldChange", colors=paletteColorBrewerRdBu,  style.name = "WikiPathways")


#Create pathway maps for all DEG pathways. Most useful 

lapply(ewp.de.wpids, function (x) {
  commandsRun(paste0('wikipathways import-as-pathway id=',x))
  loadTableData(Cortex.expr, data.key.column = "X", table.key.column = "Ensembl")
  toggleGraphicsDetails()
  setNodeColorMapping("log2FoldChange", colors=paletteColorBrewerRdBu,  style.name = "WikiPathways")
})



