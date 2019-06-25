if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8")
library(Biobase)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")
library(GEOquery)

install.packages("reshape2")
library(reshape2)

install.packages("ggplot2")
library(ggplot2)

install.packages("ggrepel")
library(ggrepel)

## load Gene microarray data
gset <- getGEO("GSE37794", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL5175", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

## Extract the expression data
edata <- exprs(gset)

## Extract the phenotype information
pdata <- pData(gset)

## Extract the features (probes) information
fdata <- fData(gset)

## Transposing the expression data such that the rows are the samples and columns are the probes
edata2 <- t(edata)

## Adding the Sample type information to be able to group all the control and DM2 samples
data.DF3 <- data.frame(Sample.Type = pdata$source_name_ch1, edata2)

## Separating the control and DM2 samples into separate data frames.
data.control <- subset(data.DF3, Sample.Type %in% "biceps brachii, CONTROL")
data.DM2 <- subset(data.DF3, Sample.Type %in% "biceps brachii,DM2")

## Removing the sample type information so as to convert the entire data frame to numerical values
data.control$Sample.Type <- NULL
data.DM2$Sample.Type <- NULL

## Calculating the mean expression values for each gene across the control and DM2 samples
data.control["Mean.Control",] <- colMeans(data.control)
data.DM2["Mean.DM2",] <- colMeans(data.DM2)

## Transposing the previous data frames so that it is easier to access the mean values as a column
data.control.2 <- as.data.frame(t(data.control))
data.DM2.2 <- as.data.frame(t(data.DM2))

## Creating a data frame for plotting which would include the individual mean values, log2 fold ratios and the p-values. 
data.plot <- data.frame(Probe.ID = rownames(data.control.2),
                        Mean.Control = data.control.2$Mean.Control,
                        Mean.DM2 = data.DM2.2$Mean.DM2)

## Log2FC calculation
data.plot$log2FC <- with(data.plot, log2(Mean.DM2/Mean.Control))

## Transpose the data.DF3 after removing the Sample type information to remove the non-numerical values, so that the probes are across the rows
data.DF4 <- data.DF3
data.DF4$Sample.Type <- NULL
data.DF4 <- as.data.frame(t(data.DF4))

## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[1:10]
  D <- x[11:20]
  p <- t.test(C,D, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
data.plot$pval <- apply(X = data.DF4, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
data.plot$log10.pval <- -log10(data.plot$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
data.plot$p.adj <- p.adjust(data.plot$pval, method = "fdr")
data.plot$log10.padj <- -log10(data.plot$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(data.plot, aes(log2FC, log10.padj, label = Probe.ID)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-0.8,0.8)) +           # Set the X-axis limits such the plot is symmetrical
  geom_point(data = subset(data.plot, log10.padj > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(data.plot, log10.padj > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(data.plot, log10.padj > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +    # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  geom_text_repel(data = subset(data.plot, log10.padj > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(DM2/CONTROL)") + ylab("-log10(FDR adjusted p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of GSE37794",
          subtitle = "With FDR Corrected p-values") +     # Plot Title
  theme_bw() +  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) # This is to center the plot title. Usually it is left-indented
ggsave(filename = "VolcanoPlot_padj_DE_Analysis_37794.png", width = 6, height = 6,
       units = "in", dpi = 300, device = "png")