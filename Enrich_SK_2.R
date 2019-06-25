install.packages("enrichR")
library(enrichR)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.8")
library(org.Hs.eg.db)

install.packages("dplyr")
library(dplyr)

## Input the list of differentially expressed genes. 
DE.genes <- c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr")

# ## To check the list of available databases (Informational purpose only. Need NOT run this line)
# dbs <- listEnrichrDbs()

## Available GO and pathway enrichment databases.
## GO_**_**_2015 databases are Gene Ontology sets for genes, ChEA is for transcription factors, KEGG is pathway enrichment database for genes and miRNAs. 
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015" , "ChEA_2016" ,"KEGG_2016")

## Enrichment Analysis
enriched <- enrichr(DE.genes, dbs)

## Write the top 10 pathways or Ontology terms for each of the databases in dbs. 
printEnrich(enriched, "EA_output.txt" , sep = "\t", columns = c(1:9))


## Visualizations Section

kgXref <- read.delim("kgXref.txt", sep = "\t", header = F, stringsAsFactors = F, na.strings = "")
names(kgXref) <- c("kgID","mRNA","spID","spDisplayID","gene.Symbol",
                   "refseq.ID","protAcc","Gene.Name","rfamAcc","tRNA.Name")

kgXref <- kgXref[,c("gene.Symbol","refseq.ID","kgID")]


## List of DE genes
DE.genes <- c("ARG2", "FHOD1", "HPGD", "VTI1B")

## List of all genes
BG.genes <- unique(kgXref$gene.Symbol)

gene.df <- bitr(DE.genes, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

universe.df <- bitr(BG.genes, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)


##### KEGG Pathway Analysis #####

KEGG.analysis <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = 'hsa',
                 universe = universe.df$ENTREZID,
                 maxGSSize = 500,
                 minGSSize = 10,
                 pvalueCutoff = 0.05)

print(KEGG.analysis)

barplot(KEGG.analysis, showCategory=20)

dotplot(KEGG.analysis, showCategory=20)




##### GO Enrichment Analysis #####

## Replace ont = “BP” with “MF” or “CC” appropriately

GO.analysis <- enrichGO(gene = gene.df$ENTREZID,
                universe = universe.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                maxGSSize = 500,
                minGSSize = 10,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable = TRUE)

print(GO.analysis)

barplot(GO.analysis, showCategory=30)

dotplot(GO.analysis, showCategory=30)

