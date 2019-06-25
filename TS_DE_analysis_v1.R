install.packages("dplyr")
library(dplyr)

fname <- "Predicted_Targets_Context_Scores.default_predictions.txt"

data <- read.table(fname, sep = "\t", header = T, na.strings = "", stringsAsFactors = FALSE)


TSdata.HS <- subset(data, Gene.Tax.ID == 9606)


kgXref <- read.delim("kgXref.txt", sep = "\t", header = F, stringsAsFactors = F, na.strings = "")
names(kgXref) <- c("kgID","mRNA","spID","spDisplayID","gene.Symbol","refseq.ID","protAcc","Gene.Name","rfamAcc","tRNA.Name")


kgXref <- kgXref[,c("gene.Symbol","refseq.ID","kgID")]

data.plot$Probe.ID <- as.numeric(sub(pattern = "X",replacement = "", x = as.character(data.plot$Probe.ID)))

data.plot2 <- merge(fdata[,c("ID","GB_LIST")], data.plot, by.x = "ID", by.y = "Probe.ID", all.x = F, all.y = F)

# Suggested by 'Chinmay Patil', https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows
S <- strsplit(data.plot2$GB_LIST, split = ",")
data.plot3 <- data.frame(ID = rep(data.plot2$ID, sapply(S, length)), RefSeq.ID = unlist(S),
                         Mean.Control = rep(data.plot2$Mean.Control, sapply(S, length)),
                         Mean.DM2 = rep(data.plot2$Mean.DM2, sapply(S, length)),
                         log2FC = rep(data.plot2$log2FC, sapply(S, length)),
                         pval = rep(data.plot2$pval, sapply(S, length)),
                         log10.pval = rep(data.plot2$log10.pval, sapply(S, length)),
                         p.adj = rep(data.plot2$p.adj, sapply(S, length)),
                         log10.padj = rep(data.plot2$log10.padj, sapply(S, length)))



data.plot4 <- unique(merge(kgXref, data.plot3, by.x = "refseq.ID", by.y = "RefSeq.ID",
                    all.x = F, all.y = T))


TSdata.HS.2 <- unique(merge(kgXref, TSdata.HS, by.x = "kgID", by.y = "Transcript.ID",
                            all.x = F, all.y = F))


data.plot4$RefSeq.Label <- substr(data.plot4$refseq.ID, start = 1, stop = 2)

data.plot5 <- subset(data.plot4, RefSeq.Label %in% c("NM","NR"))




DE.dataplot <- subset(data.plot5, abs(log2FC) >= log2(1.15) & p.adj <= 0.05)



DE.TSdata <- subset(TSdata.HS, Gene.Symbol %in% DE.dataplot$gene.Symbol)

write.table(x = DE.TSdata, file = "DE_TSdata.txt", sep = "\t", row.names = F, col.names = T)
