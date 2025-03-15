BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("DEGreport")
BiocManager::install("factoextra")
install.packages("pheatmap")
install.packages("umap")


library(DESeq2)
library(tidyverse)
library(dplyr)
library("pheatmap")


## pre-process

# load data and transpose (since we want rows as gene names)
count_data <- t(read.csv("/Users/tanlongzhi/R/plate-c/rna_ci994_241007a_renamed.csv", row.names = 1))

# create metadata
meta_data <- data.frame(sample_name = colnames(count_data))
rownames(meta_data) <- colnames(count_data)
meta_data <- meta_data %>% separate(sample_name, "_", into = c("batch", "condition", "replicate_name"))
table(meta_data$condition)

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta_data, design = ~ batch + condition)

# filter genes that are rarely expressed
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

# write DESeq2-normalized data for plotting and further analysis in MATLAB
rld <- rlog(dds, blind=FALSE)
#vsd <- vst(dds, blind=FALSE) # alternative
write.table(assay(rld), file="~/R/plate-c/rna_ci994_241007a_deseq2_matrix.csv", sep = ",", row.names=TRUE, col.names=TRUE, quote = FALSE)


## identify DEGs with Wald tests

# set significance thresholds
alpha_level=0.01
fc_threshold=0.5
#fc_threshold=0 # alternative

# differentiation: vehicle 24h vs. 0h
res_vehicle_div2_vs_div1 <- results(dds, contrast = c("condition", "DIV2-Vehicle", "DIV1-Start"), alpha=alpha_level, lfcThreshold=fc_threshold, altHypothesis="greaterAbs")
resOrdered_vehicle_div2_vs_div1 <- res_vehicle_div2_vs_div1[order(res_vehicle_div2_vs_div1$pvalue),]
summary(res_vehicle_div2_vs_div1)
write.table(cbind(gene = rownames(resOrdered_vehicle_div2_vs_div1), resOrdered_vehicle_div2_vs_div1), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div2_vs_div1.csv", sep = ",", row.names=F, col.names=TRUE, quote = FALSE)
write.table(rownames(resOrdered_vehicle_div2_vs_div1[which(resOrdered_vehicle_div2_vs_div1$log2FoldChange<fc_threshold*-1 & resOrdered_vehicle_div2_vs_div1$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div2_vs_div1.genes_down.txt", row.names=F, col.names=F, quote = FALSE)
write.table(rownames(resOrdered_vehicle_div2_vs_div1[which(resOrdered_vehicle_div2_vs_div1$log2FoldChange>fc_threshold & resOrdered_vehicle_div2_vs_div1$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div2_vs_div1.genes_up.txt", row.names=F, col.names=F, quote = FALSE)

# differentiation: vehicle 72h vs. 0h
res_vehicle_div4_vs_div1 <- results(dds, contrast = c("condition", "DIV4-Vehicle", "DIV1-Start"), alpha=alpha_level, lfcThreshold=fc_threshold, altHypothesis="greaterAbs")
resOrdered_vehicle_div4_vs_div1 <- res_vehicle_div4_vs_div1[order(res_vehicle_div4_vs_div1$pvalue),]
summary(res_vehicle_div4_vs_div1)
write.table(cbind(gene = rownames(res_vehicle_div4_vs_div1), res_vehicle_div4_vs_div1), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div4_vs_div1.csv", sep = ",", row.names=F, col.names=TRUE, quote = FALSE)
write.table(rownames(resOrdered_vehicle_div4_vs_div1[which(resOrdered_vehicle_div4_vs_div1$log2FoldChange<fc_threshold*-1 & resOrdered_vehicle_div4_vs_div1$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div4_vs_div1.genes_down.txt", row.names=F, col.names=F, quote = FALSE)
write.table(rownames(resOrdered_vehicle_div4_vs_div1[which(resOrdered_vehicle_div4_vs_div1$log2FoldChange>fc_threshold & resOrdered_vehicle_div4_vs_div1$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_vehicle_div4_vs_div1.genes_up.txt", row.names=F, col.names=F, quote = FALSE)

# HDACi: HDACi vs. vehicle at 24h
res_div2_ci994_vs_vehicle <- results(dds, contrast = c("condition", "DIV2-CI-994", "DIV2-Vehicle"), alpha=alpha_level, lfcThreshold=fc_threshold, altHypothesis="greaterAbs")
resOrdered_div2_ci994_vs_vehicle <- res_div2_ci994_vs_vehicle[order(res_div2_ci994_vs_vehicle$pvalue),]
summary(res_div2_ci994_vs_vehicle)
write.table(cbind(gene = rownames(res_div2_ci994_vs_vehicle), res_div2_ci994_vs_vehicle), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div2.csv", sep = ",", row.names=F, col.names=TRUE, quote = FALSE)
write.table(rownames(resOrdered_div2_ci994_vs_vehicle[which(resOrdered_div2_ci994_vs_vehicle$log2FoldChange<fc_threshold*-1 & resOrdered_div2_ci994_vs_vehicle$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div2.genes_down.txt", row.names=F, col.names=F, quote = FALSE)
write.table(rownames(resOrdered_div2_ci994_vs_vehicle[which(resOrdered_div2_ci994_vs_vehicle$log2FoldChange>fc_threshold & resOrdered_div2_ci994_vs_vehicle$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div2.genes_up.txt", row.names=F, col.names=F, quote = FALSE)

# HDACi: HDACi vs. vehicle at 72h
res_div4_ci994_vs_vehicle <- results(dds, contrast = c("condition", "DIV4-CI-994", "DIV4-Vehicle"), alpha=alpha_level, lfcThreshold=fc_threshold, altHypothesis="greaterAbs")
resOrdered_div4_ci994_vs_vehicle <- res_div4_ci994_vs_vehicle[order(res_div4_ci994_vs_vehicle$pvalue),]
summary(res_div4_ci994_vs_vehicle)
write.table(cbind(gene = rownames(res_div4_ci994_vs_vehicle), res_div4_ci994_vs_vehicle), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div4.csv", sep = ",", row.names=F, col.names=TRUE, quote = FALSE)
write.table(rownames(resOrdered_div4_ci994_vs_vehicle[which(resOrdered_div4_ci994_vs_vehicle$log2FoldChange<fc_threshold*-1 & resOrdered_div4_ci994_vs_vehicle$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div4.genes_down.txt", row.names=F, col.names=F, quote = FALSE)
write.table(rownames(resOrdered_div4_ci994_vs_vehicle[which(resOrdered_div4_ci994_vs_vehicle$log2FoldChange>fc_threshold & resOrdered_div4_ci994_vs_vehicle$padj < alpha_level),]), file="~/R/plate-c/rna_ci994_241007a_deseq2_log2fc_ci994_vs_vehicle_div4.genes_up.txt", row.names=F, col.names=F, quote = FALSE)
