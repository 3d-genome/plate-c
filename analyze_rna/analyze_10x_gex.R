install.packages("Matrix")
install.packages('Seurat')
install.packages('devtools')
devtools::install_github('immunogenomics/presto') # for faster marker gene analysis

library(Seurat)
library(dplyr)
library(patchwork)

## pre-process and integrate

# read 4 samples from cell ranger output
dmso_female.data <- Read10X(data.dir = "/Users/tanlongzhi/Research/dip-c/novogene/bcl-convert_sample_sheets/X202SC24072068-Z01-F008/10x/parasar_240914a_control_female_01")
dmso_male.data <- Read10X(data.dir = "/Users/tanlongzhi/Research/dip-c/novogene/bcl-convert_sample_sheets/X202SC24072068-Z01-F008/10x/parasar_240914a_control_male_01")
tsa_female.data <- Read10X(data.dir = "/Users/tanlongzhi/Research/dip-c/novogene/bcl-convert_sample_sheets/X202SC24072068-Z01-F008/10x/parasar_240914a_tsa_female_01")
tsa_male.data <- Read10X(data.dir = "/Users/tanlongzhi/Research/dip-c/novogene/bcl-convert_sample_sheets/X202SC24072068-Z01-F008/10x/parasar_240914a_tsa_male_01")

# create seurat objects
dmso_female <- CreateSeuratObject(counts = dmso_female.data, project = "dmso_female")
dmso_male <- CreateSeuratObject(counts = dmso_male.data, project = "dmso_male")
tsa_female <- CreateSeuratObject(counts = tsa_female.data, project = "tsa_female")
tsa_male <- CreateSeuratObject(counts = tsa_male.data, project = "tsa_male")

# merge into one seurat object with 4 layers
brain <- merge(dmso_female, y = c(dmso_male, tsa_female, tsa_male), add.cell.ids = c("dmso_female", "dmso_male", "tsa_female", "tsa_male"), project = "brain")

head(colnames(brain))
table(brain$orig.ident)

# analysis without integration
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
brain <- ScaleData(brain)
brain <- RunPCA(brain)
brain <- FindNeighbors(brain, dims = 1:30, reduction = "pca")
brain <- FindClusters(brain, resolution = 2, cluster.name = "unintegrated_clusters")
brain <- RunUMAP(brain, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(brain, reduction = "umap.unintegrated", group.by = c("orig.ident"))
saveRDS(brain, file = "~/R/plate-c/tsa_dmso_p6_unintegrated_240917a.RDS")

# integration
brain <- IntegrateLayers(object = brain, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
#brain <- IntegrateLayers(object = brain, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony") # alternative: harmony
brain[["RNA"]] <- JoinLayers(brain[["RNA"]])
brain <- FindNeighbors(brain, reduction = "integrated.cca", dims = 1:30)
#brain <- FindNeighbors(brain, reduction = "harmony", dims = 1:30) # alternative: harmony
brain <- FindClusters(brain, resolution = 1)
brain <- RunUMAP(brain, dims = 1:30, reduction = "integrated.cca")
#brain <- RunUMAP(brain, dims = 1:30, reduction = "harmony") # alternative: harmony
saveRDS(brain, file = "~/R/plate-c/tsa_dmso_p6_integrated_240917a.RDS")

DimPlot(brain, reduction = "umap", group.by = c("orig.ident"))
DimPlot(brain, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)

## annotate cell types

# load saved data
brain <- readRDS('~/R/plate-c/tsa_dmso_p6_integrated_240917a.RDS')
#brain_unintegrated <- readRDS("~/R/plate-c/tsa_dmso_p6_unintegrated_240917a.RDS")
#brain_harmony <- readRDS('~/R/plate-c/tsa_dmso_p6_integrated_harmony_240917a.RDS') # alternative: harmony

# visualize expression of marker genes
FeaturePlot(brain, reduction = "umap", features = c("P2ry12", "Cx3cr1", "Itgam","Mrc1","Ly86"))  # microglia
FeaturePlot(brain, reduction = "umap", features = c("Pcp2", "Itpr1"))  # purkinje
FeaturePlot(brain, reduction = "umap", features = c("Drd1", "Drd2", "Ppp1r1b", "Isl1", "Rarb"))  # msn
FeaturePlot(brain, reduction = "umap", features = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a3"))  # astro
FeaturePlot(brain, reduction = "umap", features = c("Gria1","Slc7a10","Prdm16","Gabrg3"))  # astro
FeaturePlot(brain, reduction = "umap", features = c("Abcc9","Pdgfrb", "Rgs5")) #vsmc/pericyte
FeaturePlot(brain, reduction = "umap", features = c("Flt1","Cldn5", "Pecam1", "Vwf", "Flt1", "Kdr"))  # endo
FeaturePlot(brain, reduction = "umap", features = c("Pdgfrb", "Cxcl12", "Col1a1", "Col3a1", "Col5a2", "Pdgfra"))  # vlmc
FeaturePlot(brain, reduction = "umap", features = c("Foxj1", "Dnah5", "Tmem212"))  # ependymal
FeaturePlot(brain, reduction = "umap", features = c("Kcnj13", "Ttr"))  # cp
FeaturePlot(brain, reduction = "umap", features = c("Cdh1", "Slc47a1"))  # abc
FeaturePlot(brain, reduction = "umap", features = c("Pdgfra", "Sox10", "Olig2"))  # opc
FeaturePlot(brain, reduction = "umap", features = c("Mbp", "Plp1", "Cnp", "Mag", "Sox10"))  # oligo
FeaturePlot(brain, reduction = "umap", features = c("Gabra6","Gli2","Dcc"))  # granule
FeaturePlot(brain, reduction = "umap", features = c("Tle4"))  # l6
FeaturePlot(brain, reduction = "umap", features = c("Cux2"))  # l2/3
FeaturePlot(brain, reduction = "umap", features = c("Rorb"))  # l4
FeaturePlot(brain, reduction = "umap", features = c("Nr4a2"))  # clau
FeaturePlot(brain, reduction = "umap", features = c("Prima1"))  # oec
FeaturePlot(brain, reduction = "umap", features = c("Slc17a7", "Camk2a", "Grin1", "Grin2a","Chat"))  # exc
FeaturePlot(brain, reduction = "umap", features = c("Slc17a7","Slc17a6","Gad1", "Gad2"))  # inh vs exc
FeaturePlot(brain, reduction = "umap", features = c("Ebf1","Ebf3","Zic1","Zic4"))  # granule progenitor

# merge and annotate clusters, called "cell_types_split"
brain_renamed <- RenameIdents(
  brain,
  "64" = "olfactory ensheathing cell",
  "34" = "microglia",
  "39" = "microglia",
  "27" = "endothelial cell",
  "48" = "choroid plexus epithelial cell",
  "42" = "vascular smooth muscle cell/pericyte",
  "30" = "vascular and leptomeningeal cell",
  "50" = "vascular and leptomeningeal cell",
  "57" = "arachnoid barrier cell",
  "51" = "ependymal cell",
  "6" = "oligodendrocyte progenitor",
  "55" = "oligodendrocyte progenitor",
  "32" = "oligodendrocyte progenitor",
  "28" = "oligodendrocyte",
  "31" = "oligodendrocyte",
  "8" = "astrocyte",
  "9" = "astrocyte",
  "25" = "astrocyte",
  "29" = "astrocyte",
  "12" = "astrocyte",
  "63" = "astrocyte",
  "59" = "purkinje cell",
  "7" = "medium spiny neuron",
  "17" = "medium spiny neuron",
  "43" = "medium spiny neuron",
  "10" = "cerebellar interneuron",
  "18" = "cerebellar granule cell progenitor",
  "3" = "cerebellar granule cell progenitor",
  "2" = "cerebellar granule cell",
  "5" = "cerebellar granule cell",
  "33" = "hippocampal granule cell",
  "22" = "other excitatory neuron",
  "26" = "other excitatory neuron",
  "20" = "other excitatory neuron",
  "14" = "other excitatory neuron",
  "15" = "other excitatory neuron",
  "19" = "other excitatory neuron",
  "16" = "other excitatory neuron",
  "13" = "other excitatory neuron",
  "36" = "other excitatory neuron",
  "40" = "other excitatory neuron",
  "47" = "other excitatory neuron",
  "46" = "other excitatory neuron",
  "4" = "other excitatory neuron",
  "53" = "other excitatory neuron",
  "62" = "other excitatory neuron",
  "37" = "other excitatory neuron",
  "45" = "other excitatory neuron",
  "58" = "other inhibitory neuron",
  "35" = "other inhibitory neuron",
  "52" = "other inhibitory neuron",
  "21" = "other inhibitory neuron",
  "11" = "other inhibitory neuron",
  "38" = "other inhibitory neuron",
  "1" = "unknown neuron",
  "66" = "unknown neuron",
  "60" = "unknown neuron",
  "67" = "unknown neuron",
  "65" = "unknown neuron",
  "24" = "unknown neuron",
  "49" = "unknown neuron",
  "61" = "unknown neuron",
  "0" = "unknown",
  "54" = "unknown",
  "23" = "unknown",
  "41" = "unknown",
  "56" = "unknown",
  "44" = "unknown"
)
DimPlot(brain_renamed, reduction = "umap", label = TRUE)
brain_renamed[["cell_types_split"]] <- Idents(brain_renamed)

# merge cerebellar granule cells and prognitors, called "cell_types" (the main annotation)
brain_renamed <- RenameIdents(brain_renamed, 
                              "cerebellar granule cell" = "cerebellar granule cell/progenitor",
                              "cerebellar granule cell progenitor" = "cerebellar granule cell/progenitor"
)
DimPlot(brain_renamed, reduction = "umap", label = TRUE)
brain_renamed[["cell_types"]] <- Idents(brain_renamed)

# save annotated Seurat object
saveRDS(brain_renamed, file = "~/R/plate-c/tsa_dmso_p6_integrated_240917a_annotated_250214a.RDS")
brain <- brain_renamed

## write Seurat output for plotting in MATLAB

# load saved data
brain <- readRDS('~/R/plate-c/tsa_dmso_p6_integrated_240917a_annotated_250214a.RDS')
brain_unintegrated <- readRDS("~/R/plate-c/tsa_dmso_p6_unintegrated_240917a.RDS")
#brain_harmony <- readRDS('~/R/plate-c/tsa_dmso_p6_integrated_harmony_240917a.RDS') # alternative: harmony

# write UMAP
write.table(Embeddings(brain_unintegrated, "umap.unintegrated"), "umap_unintegrated.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(Embeddings(brain, "umap"), "umap.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(as.numeric(factor(brain_unintegrated$orig.ident)), "orig.ident.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
levels(factor(brain_unintegrated$orig.ident))
DimPlot(brain_unintegrated, reduction = "umap.unintegrated", group.by = c("orig.ident"))

write.table(as.numeric(factor(brain$cell_types)), "cell_types.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.numeric(factor(brain$cell_types_split)), "cell_types_split.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
levels(factor(brain$cell_types))
levels(factor(brain$cell_types_split))

# write Harmony UMAP
#write.table(Embeddings(brain_harmony, "umap"), "umap_harmony.txt", sep = "\t", row.names = FALSE, col.names = FALSE) # alternative: harmony

# write expression vector of a gene of interest
gene_name <- "Crppa"
write.table(brain[["RNA"]]$data[gene_name,], paste0('expression_',gene_name,'.txt'),sep = "\t", row.names = FALSE, col.names = FALSE)


## make a dot plot for marker genes across cell types
DotPlot(brain, features =c("Slc1a3","Aldh1l1","Pdgfra","Sox10","Mbp","Mag","Cx3cr1","P2ry12","Dnah5","Tmem212","Kcnj13","Ttr","Frzb","Prima1","Col1a1","Col3a1","Pecam1","Kdr","Abcc9","Rgs5","Cdh1","Slc47a1","Ebf3","Boc","Gabra6","Adamts18","Rarb","Bcl11b","Pax2","Kit","Prox1","Nrp1","Nrp2","Pcp2","Calb1","Slc17a6","Slc17a7","Gad1","Gad2"),cols=c("lightgrey","red"),group.by="cell_types_split",cluster.idents=TRUE,col.min=0,col.max=2)


## prepare for DEG analysis

# load saved data
brain <- readRDS(file = "~/R/plate-c/tsa_dmso_p6_integrated_240917a_annotated_250214a.RDS")

# load list of protein-coding genes
protein_coding_genes <- readLines("/Users/tanlongzhi/R/plate-c/refdata-gex-GRCm39-2024-A.gene_names.protein_coding.txt")

# create pseudobulk object with "cell_types" annotation
brain$cell_types.orig.ident <- paste(brain$cell_types, brain$orig.ident, sep = "_")
Idents(brain) <- "cell_types.orig.ident"
pseudo_brain <- AggregateExpression(brain, return.seurat = T, group.by = c("cell_types", "orig.ident"))
table(Idents(pseudo_brain))

# create pseudobulk object with "cell_types_split" annotation (for differentiation DEGs)
brain$cell_types_split.orig.ident <- paste(brain$cell_types_split, brain$orig.ident, sep = "_")
Idents(brain) <- "cell_types_split.orig.ident"
pseudo_brain_split <- AggregateExpression(brain, return.seurat = T, group.by = c("cell_types_split", "orig.ident"))
table(Idents(pseudo_brain_split))

# create pseudobulk object with all cells
bulk_brain <- AggregateExpression(brain, return.seurat = T, group.by = "orig.ident")


## identify HDACi DEGs in all-cell aggregated bulk

# identify DEGs
degs_bulk <- FindMarkers(bulk_brain, ident.1 = c("tsa-female", "tsa-male"), ident.2 = c("dmso-female", "dmso-male"), feature=protein_coding_genes, logfc.threshold=0.5, test.use="DESeq2", min.cells.group=2, verbose = FALSE)
degs_up <- degs_bulk[which(degs_bulk$p_val_adj < 0.01 & degs_bulk$avg_log2FC > 0.5),]
degs_down <- degs_bulk[which(degs_bulk$p_val_adj < 0.01 & degs_bulk$avg_log2FC < -0.5),]
print(paste0('DEGs for full bulk: # Up: ',nrow(degs_up),'. # Down: ', nrow(degs_down),'. # Total: ',nrow(degs_up)+nrow(degs_down)))

# write DEGs to files
write.csv(degs_up, paste0('~/R/plate-c/degs_up.all.csv'))
write.csv(degs_down, paste0('~/R/plate-c/degs_down.all.csv'))

# visualize all-cell pseudobulk DEGs
genes.to.label <- c(rownames(degs_up),rownames(degs_down))
p1 <- CellScatter(bulk_brain, "dmso-female", "dmso-male", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p3 <- CellScatter(bulk_brain, "tsa-female", "tsa-male", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
p5 <- CellScatter(bulk_brain, "dmso-female", "tsa-female", highlight = genes.to.label)
p6 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE)
p7 <- CellScatter(bulk_brain, "dmso-male", "tsa-male", highlight = genes.to.label)
p8 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE)
p2 + p4 + p6 + p8

# visualize all-cell pseudobulk expression matrix for plotting in MATLAB
write.table(bulk_brain[["RNA"]]$data[protein_coding_genes,], 'pseudobulk_all.txt',sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(bulk_brain[["RNA"]]$data[c(rownames(degs_up),rownames(degs_down)),], 'pseudobulk_all_degs.txt',sep = "\t", row.names = FALSE, col.names = FALSE)


## identify HDACi DEGs in a specific cell type

# specific a cell type
deg_cell_type <- "cerebellar granule cell/progenitor"
#deg_cell_type <- "purkinje cell"
#deg_cell_type <- "medium spiny neuron"
#deg_cell_type <- "astrocyte"
#deg_cell_type <- "cerebellar interneuron"
#deg_cell_type <- "other excitatory neuron"
length(which(brain$cell_types == deg_cell_type)) # display total number of cells

# identify DEGs by single-cell (sc) and pseudobulk tests
Idents(brain) <- "cell_types.orig.ident"
degs_sc <- FindMarkers(brain, ident.1 = c(paste0(deg_cell_type,"_tsa_female"), paste0(deg_cell_type,"_tsa_male")), ident.2 = c(paste0(deg_cell_type,"_dmso_female"), paste0(deg_cell_type,"_dmso_male")), feature=protein_coding_genes, logfc.threshold=0.5, verbose = FALSE)
degs_pseudo <- FindMarkers(pseudo_brain, ident.1 = c(paste0(deg_cell_type,"_tsa-female"), paste0(deg_cell_type,"_tsa-male")), ident.2 = c(paste0(deg_cell_type,"_dmso-female"), paste0(deg_cell_type,"_dmso-male")), feature=protein_coding_genes, logfc.threshold=0.5, test.use="DESeq2", min.cells.group=2, verbose = FALSE)

# merge single-cell and pseudobulk results
names(degs_pseudo) <- paste0(names(degs_pseudo), ".bulk")
degs_pseudo$gene <- rownames(degs_pseudo)
names(degs_sc) <- paste0(names(degs_sc), ".sc")
degs_sc$gene <- rownames(degs_sc)
merge_dat <- merge(degs_sc, degs_pseudo, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# identify DEGs by intersecting single-cell and pseudobulk tests
degs_up <- merge_dat[which(merge_dat$p_val_adj.bulk < 0.01 & merge_dat$p_val_adj.sc < 0.01 & merge_dat$avg_log2FC.sc > 0.5 & merge_dat$avg_log2FC.bulk > 0.5),]
degs_down <- merge_dat[which(merge_dat$p_val_adj.bulk < 0.01 & merge_dat$p_val_adj.sc < 0.01 & merge_dat$avg_log2FC.sc < -0.5 & merge_dat$avg_log2FC.bulk < -0.5),]
print(paste0('DEGs for ',deg_cell_type,': # Up: ',nrow(degs_up),'. # Down: ', nrow(degs_down),'. # Total: ',nrow(degs_up)+nrow(degs_down)))

# write DEGs to files
write.csv(degs_up, paste0('~/R/plate-c/degs_up.',gsub(' ','_',gsub('/',' and ',deg_cell_type)),'.csv'))
write.csv(degs_down, paste0('~/R/plate-c/degs_down.',gsub(' ','_',gsub('/',' and ',deg_cell_type)),'.csv'))

# visualize DEGs
genes.to.label <- c(degs_up$gene,degs_down$gene)
p1 <- CellScatter(pseudo_brain, c(paste0(deg_cell_type,"_dmso-female")), c(paste0(deg_cell_type,"_dmso-male")), highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p3 <- CellScatter(pseudo_brain, c(paste0(deg_cell_type,"_tsa-female")), c(paste0(deg_cell_type,"_tsa-male")), highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
p5 <- CellScatter(pseudo_brain, c(paste0(deg_cell_type,"_dmso-female")), c(paste0(deg_cell_type,"_tsa-female")), highlight = genes.to.label)
p6 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE)
p7 <- CellScatter(pseudo_brain, c(paste0(deg_cell_type,"_dmso-male")), c(paste0(deg_cell_type,"_tsa-male")), highlight = genes.to.label)
p8 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE)
p2 + p4 + p6 + p8


## identify differentiation DEGs in control cerebellar granule cells

# identify DEGs by single-cell (sc) and pseudobulk tests
Idents(brain) <- "cell_types_split.orig.ident"
degs_sc <- FindMarkers(brain, ident.1 = c("cerebellar granule cell_dmso_female", "cerebellar granule cell_dmso_male"), ident.2 = c("cerebellar granule cell progenitor_dmso_female", "cerebellar granule cell progenitor_dmso_male"), feature=protein_coding_genes, logfc.threshold=0.5, verbose = FALSE)
degs_pseudo <- FindMarkers(pseudo_brain_split, ident.1 = c("cerebellar granule cell_dmso-female", "cerebellar granule cell_dmso-male"), ident.2 = c("cerebellar granule cell progenitor_dmso-female", "cerebellar granule cell progenitor_dmso-male"), feature=protein_coding_genes, logfc.threshold=0.5, test.use="DESeq2", min.cells.group=2, verbose = FALSE)

# merge single-cell and pseudobulk results
names(degs_pseudo) <- paste0(names(degs_pseudo), ".bulk")
degs_pseudo$gene <- rownames(degs_pseudo)
names(degs_sc) <- paste0(names(degs_sc), ".sc")
degs_sc$gene <- rownames(degs_sc)
merge_dat <- merge(degs_sc, degs_pseudo, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# identify DEGs by intersecting single-cell and pseudobulk tests
degs_up <- merge_dat[which(merge_dat$p_val_adj.bulk < 0.01 & merge_dat$p_val_adj.sc < 0.01 & merge_dat$avg_log2FC.sc > 0.5 & merge_dat$avg_log2FC.bulk > 0.5),]
degs_down <- merge_dat[which(merge_dat$p_val_adj.bulk < 0.01 & merge_dat$p_val_adj.sc < 0.01 & merge_dat$avg_log2FC.sc < -0.5 & merge_dat$avg_log2FC.bulk < -0.5),]
print(paste0('DEGs for differentiation: # Up: ',nrow(degs_up),'. # Down: ', nrow(degs_down),'. # Total: ',nrow(degs_up)+nrow(degs_down)))

# write DEGs to files
write.csv(degs_up, paste0('~/R/plate-c/degs_up.differentiation.csv'))
write.csv(degs_down, paste0('~/R/plate-c/degs_down.differentiation.csv'))


## compare HDACi DEGs between two cell types

# specify cell type 1
deg_cell_type_1 <- "cerebellar granule cell/progenitor"
length(which(brain$cell_types == deg_cell_type_1))

# specify cell type 2
deg_cell_type_2 <- "medium spiny neuron"
#deg_cell_type_2 <- "purkinje cell"
#deg_cell_type_2 <- "astrocyte"
#deg_cell_type_2 <- "cerebellar interneuron"
length(which(brain$cell_types == deg_cell_type_2))

# identify pseudobulk DEGs
degs_pseudo_1 <- FindMarkers(pseudo_brain, ident.1 = c(paste0(deg_cell_type_1,"_tsa-female"), paste0(deg_cell_type_1,"_tsa-male")), ident.2 = c(paste0(deg_cell_type_1,"_dmso-female"), paste0(deg_cell_type_1,"_dmso-male")), feature=protein_coding_genes, logfc.threshold=0.5, test.use="DESeq2", min.cells.group=2, verbose = FALSE)
degs_pseudo_2 <- FindMarkers(pseudo_brain, ident.1 = c(paste0(deg_cell_type_2,"_tsa-female"), paste0(deg_cell_type_2,"_tsa-male")), ident.2 = c(paste0(deg_cell_type_2,"_dmso-female"), paste0(deg_cell_type_2,"_dmso-male")), feature=protein_coding_genes, logfc.threshold=0.5, test.use="DESeq2", min.cells.group=2, verbose = FALSE)

# merge results from both cell types
names(degs_pseudo_1) <- paste0(names(degs_pseudo_1), ".", deg_cell_type_1)
degs_pseudo_1$gene <- rownames(degs_pseudo_1)
names(degs_pseudo_2) <- paste0(names(degs_pseudo_2), ".", deg_cell_type_2)
degs_pseudo_2$gene <- rownames(degs_pseudo_2)
merge_dat <- merge(degs_pseudo_1, degs_pseudo_2, by = "gene")

# write results to files
write.csv(merge_dat, paste0('~/R/plate-c/degs.between_cell_types.',gsub(' ','_',gsub('/',' and ',deg_cell_type_1)),'.vs.',gsub(' ','_',gsub('/',' and ',deg_cell_type_2)),'.csv'))
