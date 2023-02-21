library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

count <- Read10X("Data/PM_PS_0001_T_A1/outs/filtered_feature_bc_matrix/")
count[1:4, 1:4]

myObj <- CreateSeuratObject(counts = count, project = "CRC")
myObj@meta.data$Library <- "PM-PS-0001-T-A1"
head(myObj@meta.data, n = 3)

# do not run
#annotation <- read.delim('table/PM_PS_0001_T_A1_anno.txt', row.names = 1)
#head(annotation, n = 3)
#
#raw_count <- Read10X("Data/PM_PS_0001_T_A1/outs/raw_feature_bc_matrix/")
#raw_count <- raw_count[, which(colnames(raw_count) %in% rownames(annotation))]
#
#myObj_w_anno <- CreateSeuratObject(counts = raw_count, project = "CRC")
#myObj_w_anno <- AddMetaData(myObj_w_anno, metadata = annotation)
#head(myObj_w_anno@meta.data, n = 3)
#remove(raw_count)

plt1 <- ggplot(myObj@meta.data, aes(Library, nCount_RNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
plt2 <- ggplot(myObj@meta.data, aes(Library, nCount_RNA)) + 
  geom_jitter(size = 0.25) + 
  theme_bw() + ylim(0, 2000) + 
  theme(panel.grid = element_blank())

plt3 <- ggplot(myObj@meta.data, aes(Library, nFeature_RNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
plt4 <- ggplot(myObj@meta.data, aes(Library, nFeature_RNA)) + 
  geom_jitter(size = 0.25) + 
  geom_hline(yintercept = 200, col = "red2") +
  theme_bw() + ylim(0, 1000) + 
  theme(panel.grid = element_blank())

CombinePlots(plots = list(plt1, plt2, plt3, plt4), ncol = 2)

dir.create("stats")
ggsave("stats/umi_genecounts.pdf", units = "cm", width = 14, height = 14)


### filter low quality cells
myObj <- subset(myObj, subset = nFeature_RNA >= 200)
nrow(myObj@meta.data)


umithre <- as.integer(mean(myObj@meta.data$nCount_RNA) + 2*sd(myObj@meta.data$nCount_RNA)); umithre
myObj <- subset(myObj, nCount_RNA < umithre)
nrow(myObj@meta.data)


myObj[["percent.mt"]] <- PercentageFeatureSet(myObj, pattern = "^MT-")

ggplot(myObj@meta.data, aes(Library, percent.mt)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  geom_hline(yintercept = 20, col = "red2") +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave("stats/percent_mt.pdf", units = "cm", width = 5, height = 6)

myObj <- subset(myObj, percent.mt < 20)
nrow(myObj@meta.data)

plt1 <- ggplot(myObj@meta.data, aes(Library, nCount_RNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
plt2 <- ggplot(myObj@meta.data, aes(Library, nCount_RNA)) + 
  geom_jitter(size = 0.25) + 
  theme_bw() + ylim(0, 2000) + 
  theme(panel.grid = element_blank())

plt3 <- ggplot(myObj@meta.data, aes(Library, nFeature_RNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
plt4 <- ggplot(myObj@meta.data, aes(Library, nFeature_RNA)) + 
  geom_jitter(size = 0.25) + 
  geom_hline(yintercept = 200, col = 'red2') +
  theme_bw() + ylim(0, 1000) + 
  theme(panel.grid = element_blank())

CombinePlots(plots = list(plt1, plt2, plt3, plt4), ncol = 2)
ggsave("stats/umi_genecounts.postQc.pdf", units = "cm", width = 14, height = 14)

plt1 <- FeatureScatter(myObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plt2 <- FeatureScatter(myObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plt1, plt2))
ggsave("stats/umi_genecounts_corr.pdf", units = "cm", width = 20, height = 8)


### Normalization
myObj <- NormalizeData(myObj, normalization.method = "LogNormalize", scale.factor = 10000)

myObj <- FindVariableFeatures(myObj, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(myObj), 10)
plot1 <- VariableFeaturePlot(myObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot2
ggsave("stats/hvgs.pdf", units = "cm", width = 16, height = 8)

dir.create("tmp")

myObj <- ScaleData(myObj, vars.to.regress = c('nCount_RNA', 'percent.mt'))
myObj <- RunPCA(myObj, features = VariableFeatures(myObj), npcs = 60)

# will take a lot of time
myObj <- JackStraw(myObj, num.replicate = 100, dims = 60)
### or parallelization using Future (NormalizeData, ScaleData, JackStraw, FindMarkers, FindIntegrationAnchors, FindClusters)
#library(future)
#plan("multisession", workers = 4)
#myObj <- JackStraw(myObj, num.replicate = 100, dims = 60)

myObj <- ScoreJackStraw(myObj, dims = 1:60)
saveRDS(myObj, 'tmp/myObj.Rds')

JackStrawPlot(myObj, dims = 1:60)
ggsave("stats/jackstraw_pca.pdf", units = "cm", width = 30, height = 12)

pvals <- data.frame(myObj@reductions$pca@jackstraw$overall.p.values)
pcs_use <- pvals[pvals$Score > 0.001, 'PC'][1] -1
print (paste(c('Significant PCs:', pcs_use), collapse = ' ')) # 56
# or you can set arbitrary number of PCs; 30 PCs


myObj <- RunUMAP(myObj, dims = 1:pcs_use, reduction = "pca", n.components = 2, min.dist = 0.4, seed.use = 123) # min.dist = 0.2 ~ 0.4
saveRDS(myObj, 'tmp/myObj.Rds')

dir.create("plots")
FeaturePlot(myObj, features = "nFeature_RNA", cols = c("grey", "red"))
ggsave("plots/umap_genecount.pdf", units = "cm", width = 10, height = 8)


genes <- c("EPCAM", "KRT18", "PTPRC") # epithelial, cancer, immune cells
FeaturePlot(myObj, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave("plots/umap_markers_epimmn.pdf", units = "cm", width = 18, height = 18)

genes <- c("CD3E", "CD8A", "TNFRSF4", "CD79A", "IGHA1") # T/B cells
FeaturePlot(myObj, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave("plots/umap_markers_adaptive.pdf", units = "cm", width = 18, height = 24)

genes <- c("C1QA", "HLA-DPB1", "CPA3", "TPSAB1") # myeloid, cDC, mast cell
FeaturePlot(myObj, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave("plots/umap_markers_innate.pdf", units = "cm", width = 18, height = 18)

genes <- c("COL1A1", "VWF") # fibroblast, endothelial cell
FeaturePlot(myObj, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave("plots/umap_markers_stromal.pdf", units = "cm", width = 18, height = 8)


### Clustering
myObj <- FindNeighbors(myObj, dims = 1:pcs_use)
myObj <- FindClusters(myObj, resolution = 0.4)
saveRDS(myObj, 'tmp/myObj.Rds')

DimPlot(myObj, label = T, reduction = 'umap')
ggsave("plots/umap_res0.4.pdf", units = "cm", width = 12, height = 10)

dir.create("degs")
markers <- FindAllMarkers(myObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = "MAST")
write.table(markers, "degs/markers.res_0.3.txt", sep = "\t", quote = F, col.names = T, row.names = F)

head(markers, n = 3)
avg_logFC <- 1 # or 2
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

dehm <- DoHeatmap(object = myObj, features = top10$gene, size = 3, raster = T, draw.lines = F, label = F)
ggsave('degs/markers.res_0.4.pdf', units = 'cm', width = 40, height = 30)


### Cell annotation
clusternames <- c("Epithelial cell", "T cell", "CD8 T cell", "T cell", "Tregs",
                  "Myeloid", "Fibroblast", "T cell", "T cell", "Plasma cell",
                  "B cell", "Epithelial cell", "Epithelial cell", "Endothelial cell", "pDC",
                  "Mast cell")

myObj@meta.data$Celltype <- mapvalues(myObj@meta.data$RNA_snn_res.0.4, from = c(0:15), to = clusternames)
myObj@meta.data$Celltype <- factor(myObj@meta.data$Celltype, 
                                   levels = c("Epithelial cell", "Myeloid", "pDC", "Mast cell", "B cell", 
                                              "T cell", "CD8 T cell", "Tregs", "Fibroblast", "Plasma cell", "Endothelial cell"))
summary(myObj@meta.data$Celltype)
Idents(myObj) <- "Celltype"
saveRDS(myObj, 'tmp/myObj.Rds')


DimPlot(myObj, label = T, reduction = 'umap')
ggsave("plots/umap_celltype.pdf", units = "cm", width = 14, height = 10)

markers_ct <- FindAllMarkers(myObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = "MAST")
write.table(markers_ct, "degs/markers.celltype.txt", sep = "\t", quote = F, col.names = T, row.names = F)

top10 <- markers_ct %>% group_by(cluster) %>% top_n(10, avg_logFC)

dehm <- DoHeatmap(object = myObj, features = top10$gene, size = 3, raster = T, draw.lines = F, label = F)
ggsave('degs/markers.celltype.pdf', units = 'cm', width = 40, height = 30)

