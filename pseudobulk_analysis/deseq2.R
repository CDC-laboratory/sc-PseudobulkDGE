#! /usr/bin/env Rscript
library("optparse")
 
option_list = list(
  make_option(c("-g", "--group"), type="character", 
              help="comparison group folder", metavar="character"),
  make_option(c("-r", "--contrast"), type="character", default=NULL, 
              help="contrast folder", metavar="character"),
  make_option(c("-c", "--control"), type="character", 
              help="control variable", metavar="character"),
  make_option(c("-t", "--treatment"), type="character", 
              help="treatment variable", metavar="character"),
  make_option(c("-m", "--mincount"), type="character", default=10, 
              help="min. count for gene filtering", metavar="integer"),
  make_option(c("-v", "--covariate"), type="character", default=NULL, 
              help="covariates (comma separated)", metavar="character")

);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)
print(length(opt))
# Checking if there are three arguments provided
if (length(opt)< 3) {
  stop("At least three arguments must be supplied", call.=FALSE)
}

## load libraries
suppressMessages(library(DESeq2))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(pheatmap))

convert_numeric_to_integer <- function(df) {
  # Save row and column names
  rownames_df <- rownames(df)
  colnames_df <- colnames(df)
  
  # Convert numeric columns to integers
  df_int <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) as.integer(round(col)) else col
  }))
  
  # Restore row and column names
  rownames(df_int) <- rownames_df
  colnames(df_int) <- colnames_df
  
  return(df_int)
}


fname <- paste(opt$group, opt$contrast, sep="/")
control <- opt$control
treatment <- opt$treatment
print(fname)

m <- read.csv(paste(fname, "conditions.csv", sep="/"), sep="\t", row.names = 1)
c <- read.csv(paste(fname, "counts.csv", sep="/"), sep="\t", row.names = 1)
cf <- read.csv(paste(fname, "countsabove1000.csv", sep="/"), sep="\t", row.names = 1)
print(row.names(c))

## conver to int
#c <- convert_numeric_to_integer(c)
#cf <- convert_numeric_to_integer(cf)

# Checking if index in colData (m) and in counts.csv are in the same order
if (any(row.names(c)!=row.names(m))) {
        print("Changing order of rows...")
        c <- c[row.names(m),]
        cf <- cf[row.names(m),]
}
print(row.names(m))
print(row.names(c))
print(dim(c))
## all genes

val <- ""
print(opt$covariate)

if (!is.null(opt$covariate) && length(opt$covariate) > 0)  {
        val <- paste0(opt$covariate.strreplace(",", "+"), '+')
}
print(paste0("~", val, "condition"))
dds <- DESeqDataSetFromMatrix(countData = t(c), colData=m, 
                                design=formula(paste0("~", val, "condition")))

## before running DESeq2, we will add the order of categories given control-condition/case variable
dds$condition <- factor(dds$condition, levels = c(control, treatment))

## filter low-expressed genes
countsDDs <- counts(dds)
initrow <- nrow(countsDDs)
keep <- rowSums(countsDDs) >= opt$mincount
print(length(keep))
dds <- dds[keep,]

print(paste("n. of genes filtered by min", opt$mincount, "counts:", 
        initrow - nrow(counts(dds)), sep=" "))

# Estimate size factors with poscounts to overcome log with zero values
##Uses positive counts only to compute size factors.
#dds <- estimateSizeFactors(dds, type = "poscounts")

## run statistical test and export results
dds <- DESeq(dds)
res <- results(dds)

## order by pvalue and export counts (norm and unnormalized)
resOrdered <- res[order(res$pvalue),]
normc <- counts(dds, normalized = TRUE)
con <- counts(dds, normalized = FALSE) ## the default is false

write.csv(normc, paste(fname, "DGE_normcounts.csv", sep="/"))
write.csv(con, paste(fname, "DGE_counts.csv", sep="/"))
write.csv(resOrdered, paste(fname, "res.csv", sep="/"))

res05 <- results(dds, alpha=0.05)
res05Ordered <- res05[order(res05$pvalue),]
write.csv(subset(res05Ordered, padj < 0.05), paste(fname, "res05.csv", sep="/"))
write.csv(subset(resOrdered, padj < 0.1), paste(fname, "res1.csv", sep="/"))
write.csv(subset(res05Ordered, padj < 0.05 & log2FoldChange > 2), paste(fname, "res05-UP-LFC2.csv", sep="/"))
write.csv(subset(res05Ordered, padj < 0.05 & log2FoldChange < -2), paste(fname, "res05-DOWN-LFC2.csv", sep="/"))
res05OrderedT <- as.data.frame(res05Ordered)

## quality plots
# Plot dispersion estimates
png(file=paste(fname, "dispersion-estimates.png", sep="/"))
plotDispEsts(dds)
dev.off()

### different transformations
vsd <- vst(dds, blind=TRUE) ## recommended for large sample numbers
#rld <- rlog(dds, blind=TRUE)

## plot PCA
###  vsd normalisation method
#### condition
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png(file=paste(fname, "PCA.png", sep="/"),
                width = 1200, height = 900,
                res=200)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

#### sample_ID (to add cool palette here)
pcaData <- plotPCA(vsd, intgroup=c("sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png(file=paste(fname, "PCA-replicates.png", sep="/"),
                width = 1200, height = 900, res=200)
ggplot(pcaData, aes(PC1, PC2, color=sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()


## sample dist heatmap
library("RColorBrewer")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## select 10 samples from each group
png(file=paste(fname, "10samplesdist_heatmap.png", sep="/"))

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
         #annotation = m[, c("condition"), drop=F])
dev.off()

## count plots
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)

png(file=paste(fname, "counts_siggene.png", sep="/"))
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
dev.off()

## MA plots
png(file=paste(fname, "rawMA.png", sep="/"))
plotMA(res, ylim=c(-2,2))
dev.off()
png(file=paste(fname, "shrunkMA.png", sep="/"))
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

## make vulcano plot
## lfcShrink before?
png(file=paste(fname, "volcano.png", sep="/"),
      width = 1200, height = 900, res=100)
EnhancedVolcano(res, 
                lab = rownames(res), x = 'log2FoldChange', 
                title = 'DGE',
               subtitle = NULL,
                y = 'padj', pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 6.0,
               col=c('#353236', 'grey','#68ba7b', '#409ac9'),
                colAlpha = 0.4,
                legendPosition = "right",
                legendLabSize = 18,
                xlim = c(-10,10),
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                'p-value & Log (base 2) FC'))
dev.off()

## plot norm expression of 20 most significant genes
### res05Ordered and normc
## Order results by padj values
top20_sig_genes <- res05OrderedT %>%
        rownames_to_column(var="gene") %>%
        dplyr::pull(gene) %>%
        head(n=20)

top20_sig_norm <- data.frame(normc) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes)

## we need to map our sample names from the metadata table as the names 
## given by deseq2 are not equal
gathered_top20_sig <- top20_sig_norm %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], 
        key = "samplename", value = "normalized_counts")
print(head(gathered_top20_sig))
m$sample_id <- row.names(m)
gathered_top20_sig <- inner_join(m[, c("sample_id", "condition" )], gathered_top20_sig, 
        by = c("sample_id" = "samplename"))
print(head(gathered_top20_sig))

## plot using ggplot2
png(file=paste(fname, "20sig-genes.png", sep="/"))
ggplot(gathered_top20_sig) +
        geom_point(aes(x = gene, 
                       y = normalized_counts, 
                       color = condition), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

## also the heatmap
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")
png(file=paste(fname, "20sig-genes-heatmap.png", sep="/"))
# Run pheatmap using the metadata data frame for the annotation
pheatmap(top20_sig_norm[ , 2:length(colnames(top20_sig_norm))], 
    color = heat_colors, 
    cluster_rows = F, 
    show_rownames = T,
    #annotation = m$condition,
    border_color = NA, 
    fontsize = 10, 
    scale = "row", 
    fontsize_row = 10, 
    height = 20)   
dev.off()

## plot norm counts violin
gathered_normc <- data.frame(normc) %>%
        rownames_to_column(var = "gene") %>%
        gather(colnames(normc)[2:length(colnames(normc))], 
        key = "samplename", value = "normalized_counts")
print('gathered_norm')

gathered_normc <- inner_join(m[, c("sample_id", "condition" )], gathered_normc, 
        by = c("sample_id" = "samplename"))
print('gathered_norm')

## plot using ggplot2
png(file=paste(fname, "violin-normc.png", sep="/"))
ggplot(gathered_normc) +
        geom_violin(aes(x = sample_id, 
                       y = normalized_counts, 
                       color = condition), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Sample") +
        ylab("log10 Normalized Counts") +
        ggtitle("Normalised counts distribution") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
dev.off()