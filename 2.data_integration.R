library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

### Data integration & reference-based annotation
# load seurat object
myObj <- readRDS("tmp/myObj.Rds")

# new data
count_new <- Read10X("Data/PM_PS_0022_T_A1/outs/filtered_feature_bc_matrix/")
count_new[1:4, 1:4]

newObj <- CreateSeuratObject(counts = count_new, project = "CRC")
newObj@meta.data$Library <- "PM-PS-0022-T-A1"
head(newObj@meta.data, n = 3)

newObj <- subset(newObj, subset = nFeature_RNA >= 200)
nrow(newObj@meta.data) # 3517

umithre <- as.integer(mean(newObj@meta.data$nCount_RNA) + 2*sd(newObj@meta.data$nCount_RNA)); umithre
newObj <- subset(newObj, nCount_RNA < umithre)
nrow(newObj@meta.data) # 3385

newObj[["percent.mt"]] <- PercentageFeatureSet(newObj, pattern = "^MT-")
newObj <- subset(newObj, percent.mt < 20)
nrow(newObj@meta.data) # 2315

new_plt1 <- ggplot(newObj@meta.data, aes(Library, nCount_RNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
new_plt2 <- ggplot(newObj@meta.data, aes(Library, nCount_RNA)) + 
  geom_jitter(size = 0.25) + 
  theme_bw() + ylim(0, 2000) + 
  theme(panel.grid = element_blank())

new_plt3 <- ggplot(newObj@meta.data, aes(Library, nFeature_RNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
new_plt4 <- ggplot(newObj@meta.data, aes(Library, nFeature_RNA)) + 
  geom_jitter(size = 0.25) + 
  geom_hline(yintercept = 200, col = 'red2') +
  theme_bw() + ylim(0, 1000) + 
  theme(panel.grid = element_blank())

CombinePlots(plots = list(new_plt1, new_plt2, new_plt3, new_plt4), ncol = 2)
ggsave("stats/newdata.umi_genecounts.postQc.pdf", units = "cm", width = 14, height = 14)

newObj <- NormalizeData(newObj, normalization.method = "LogNormalize", scale.factor = 10000)
newObj <- FindVariableFeatures(newObj, selection.method = "vst", nfeatures = 2000)

myObj@meta.data$seurat_clusters <- NULL
myObj@meta.data$RNA_snn_res.0.4 <- NULL
head(myObj@meta.data, n=3)
head(newObj@meta.data, n=3)


### Label transfer
transfer.anchors <- FindTransferAnchors(reference = myObj, query = newObj, dims = 1:30)

predictions_type <- TransferData(anchorset = transfer.anchors, refdata = myObj@meta.data$Celltype, dims = 1:30)
saveRDS(predictions_type, 'tmp/predictions_type.Rds')
head(predictions_type, n=3)
nrow(subset(predictions_type, prediction.score.max >= 0.5)) # 2314


newObj_pred <- subset(newObj, cells = rownames(predictions_type))
identical(rownames(newObj_pred@meta.data), rownames(predictions_type))

newObj_pred@meta.data$Celltype <- predictions_type$predicted.id
head(newObj_pred@meta.data, n = 3)

summary(as.factor(newObj_pred@meta.data$Celltype))
saveRDS(newObj_pred, "tmp/newObj_pred.Rds")


### Merge objects
mergedObj <- merge(x = myObj, y = newObj_pred)
head(mergedObj@meta.data, n=3)

mergedObj@meta.data$Library <- factor(mergedObj@meta.data$Library, levels = c("PM-PS-0001-T-A1", "PM-PS-0022-T-A1"))
mergedObj@meta.data$Celltype <- factor(mergedObj@meta.data$Celltype, 
                                       levels = c("Epithelial cell", "Myeloid", "pDC", "Mast cell", "B cell", 
                                                  "T cell", "CD8 T cell", "Tregs", "Fibroblast", "Plasma cell", "Endothelial cell"))
Idents(mergedObj) <- "Celltype"
saveRDS(mergedObj, 'tmp/mergedObj.Rds')


mergedObj <- NormalizeData(object = mergedObj, normalization.method = "LogNormalize", scale.factor = 10000)
mergedObj <- FindVariableFeatures(object = mergedObj, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = mergedObj), 10)
plot1 <- VariableFeaturePlot(object = mergedObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot2
ggsave("stats/hvgs.merged.pdf", units = 'cm', width = 17, height = 10)

mergedObj <- ScaleData(object = mergedObj, vars.to.regress = c('nCount_RNA', 'percent.mt'))
mergedObj <- RunPCA(object = mergedObj, features = VariableFeatures(mergedObj), npcs = 70)
mergedObj <- JackStraw(object = mergedObj, num.replicate = 100, dims = 70)
mergedObj <- ScoreJackStraw(object = mergedObj, dims = 1:70)
saveRDS(mergedObj, 'tmp/mergedObj.Rds')

JackStrawPlot(mergedObj, dims = 1:70) 
ggsave("stats/jackstraw_pca.merged.pdf", units = 'cm', width = 40, height = 12)

pvals <- data.frame(mergedObj@reductions$pca@jackstraw$overall.p.values)
pcs_use <- pvals[pvals$Score > 0.001, 'PC'][1] -1
print (paste(c('Significant PCs:', pcs_use), collapse = ' ')) # 58


### Data integration - Harmony
library(harmony)

mergedObj <- RunHarmony(mergedObj, group.by.vars = "Library")

mergedObj <- RunUMAP(mergedObj, dims = 1:pcs_use, reduction = "harmony", n.components = 2, min.dist = 0.4, seed.use = 1234)
mergedObj <- RunTSNE(mergedObj, dims = 1:pcs_use, reduction = "harmony", dim.embed = 2, seed.use = 4321)
saveRDS(mergedObj, 'tmp/mergedObj.Rds')

DimPlot(mergedObj, reduction = 'umap', pt.size = .5, label = T, label.size = 2)
ggsave("plots/umap_celltype.merged.pdf", units = "cm", width = 15, height = 10)

DimPlot(mergedObj, group.by = "Library", reduction = "umap")
ggsave("plots/umap_library.merged.pdf", units = "cm", width = 16, height = 10)


DimPlot(mergedObj, reduction = "tsne", pt.size = .5, label = T, label.size = 2)
ggsave("plots/tsne_celltype.merged.pdf", units = "cm", width = 15, height = 10)

DimPlot(mergedObj, group.by = "Library", reduction = "tsne")
ggsave("plots/tsne_library.merged.pdf", units = "cm", width = 16, height = 10)


markers_ct_m <- FindAllMarkers(mergedObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = "MAST")
write.table(markers_ct_m, "degs/markers.celltype.merged.txt", sep = "\t", quote = F, col.names = T, row.names = F)

top10 <- markers_ct_m %>% group_by(cluster) %>% top_n(10, avg_logFC)

dehm <- DoHeatmap(object = mergedObj, features = top10$gene, size = 3, raster = T, draw.lines = F, label = F)
ggsave('degs/markers.celltype.merged.pdf', units = 'cm', width = 40, height = 30)

sessionInfo()

