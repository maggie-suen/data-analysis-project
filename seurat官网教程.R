####10X单细胞数据分析####

####初筛+数据初处理####
library(dplyr)
library(Seurat)
library(patchwork)
#load the PBMC data 一次提取barcodes,genes,matrix三个文件
pbmc.data <- Read10X(data.dir = "D:/1博后工作/生信学习/seurat官网教程/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#用于计算每个细胞中线粒体基因的百分比。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
# 可视化QC指标，并使用这些指标来筛选细胞
#我们筛选具有超过 2500 个或小于 200 个独特特征计数的细胞
#我们筛选线粒体计数 > 5% 的细胞
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#FeatureScatter函数
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#数据归一化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#或者这句话，同样可以归一化
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#应用线性映射(缩放), 这是PCA等降维技术之前的标准预处理步骤。ScaleData()函数
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


####PCA降维####
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways可视化定义 PCA 的单元格和特征
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#来可视化定义 PCA 的单元格和特征，包括 VizDimReduction ()、DimPlot () 和 DimHeatmap ()
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#或
DimPlot(pbmc, reduction = "pca") + NoLegend()
#或
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap中dims参数，看是否分成差异明显的群，来决定用哪些PCA主成分
#DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
#确定数据集的维度



####聚类####
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


####非线性降维####
#UMAP/tSNE两种降维方式
#UMAP降维
pbmc <- RunUMAP(pbmc, dims = 1:10)
#可以通过设置label=TRUE或使用LabelClusters函数来帮助标记单个细胞簇
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",label=TRUE,label.size = 5)
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

#tSNE降维
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
# 显示在聚类标签
DimPlot(pbmc, reduction = "tsne", label = TRUE,label.size=5)

####寻找聚类marker####
# find all markers of cluster 2
# FindAllMarkers会识别每个分簇的markers，可以测试该分簇与其他分簇或全部细胞的差异
# min.pct参数，在两组细胞中的任何一组中检测到的最小百分
# thresh.test参数，在两组细胞间以一定数量的差异表达（平均）
# max.cells.per.ident参数，通过降低每个类的采样值，提高计算速度
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2,min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3),min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE,min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# the ROC test returns the ‘classification power’ for any individual er (ranging from 0 - random, to 1 - perfect)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)




####绘图####
#VLnPlot函数,完整代码
VlnPlot(pbmc, features = c("MS4A1", "CD79A"),#基因名
        layer = "counts", #指定使用哪个数据层进行绘图
        log = TRUE, #对Y轴（表达量）进行对数转换
        pt.size = 0, #点的大小
        ncol = 2,  #分两列显示
        split.by = NULL)  #不按条件拆分


#FeaturePlot函数
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

#RidgePlot函数
RidgePlot(pbmc, features = c("MS4A1", "CD79A"))

top10 <- pbmc.ers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#DoHeatmap()为给定的单元格和特征生成表达式热图。
pbmc.markers%>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#识别细胞类型
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
#存储聚类标记图
ggsave(filename = "D:/1博后工作/生信学习/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)            
