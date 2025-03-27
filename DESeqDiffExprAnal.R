# Differential expression analysis
# chapter 5 in the url
# browseURL("http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html")

#load our DESeq2 data sets from before
ddslist<-readRDS(file = "dds.rds")
dds<-ddslist$`source_name_ch1`
  
#5.1 Running the differential expression pipeline
dds <- DESeq(dds)
#5.2 Building the results table
res <- results(dds)
res


#dds@colData$source_name_ch1
# at first, we compare IL_17A treated with control group 
res <- results(dds, contrast=c("source_name_ch1","IL_17A","CTRL"))
mcols(res, use.names = TRUE) %>% as_tibble()
summary(res)

#recall p-values and adjusted p-values: https://www.youtube.com/watch?v=K8LQSvtjcEo
res.05 <- results(dds, alpha = 0.001)
table(res.05$padj < 0.001)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.001)
#5.3 Other comparisons
#dds$source_name_ch1
#results(dds, contrast = c("characteristics_ch1.1", "timepoint..DIV.0", "timepoint..DIV.4"))

#5.4 Multiple testing
# how many genes are differentially expressed based on p-value?
#sum(res$pvalue < 0.001, na.rm=TRUE)
#sum(!is.na(res$pvalue))

# we use 0.001 threshold, based on the paper indications

sum(res$padj < 0.001, na.rm=TRUE)
# We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation:
resSig <- subset(res, padj < 0.001)
#resSig2 <- subset(resSig, abs(log2FoldChange>1))

count(resSig$log2FoldChange > 0)
count(resSig$log2FoldChange < 0)
head(resSig[ order(resSig$padj), ])
# …and with the strongest up-regulation:
head(resSig[ order(resSig$padj, decreasing = TRUE), ])


#6 Plotting results
#6.1Counts plot

#topGene <- rownames(res)[which.min(res$padj)]
topGene <- "ENSG00000162747"
plotCounts(dds, gene = topGene, intgroup=c("source_name_ch1"))
?plotCounts

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("source_name_ch1","characteristics_ch1.1"),
                         returnData = TRUE)

ggplot(geneCounts, aes(x = source_name_ch1, y = count, color = characteristics_ch1.1, group = characteristics_ch1.1)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


#MA plot

# BiocManager::install("apeglm")
# library("apeglm")

resultsNames(dds)
res <- lfcShrink(dds, coef="source_name_ch1_IL_17A_vs_CTRL", type="apeglm")
DESeq2::plotMA(res, ylim = c(-5, 5))

# gene clustering

# library("genefilter")
vsd <- vst(dds, blind = FALSE)
vsd$source_name_ch1
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
mat  <- assay(vsd)[ topVarGenes, ]
mat1 <- mat[,c(1:6,7:9,13:15)]
mat1  <- mat1 - rowMeans(mat1)
anno <- as.data.frame(colData(vsd)[,c("characteristics_ch1.1","source_name_ch1")])
pheatmap(mat = mat1, fontsize_row = 5)
colData(vsd)

