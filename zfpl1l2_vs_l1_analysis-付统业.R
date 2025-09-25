setwd("/media/dell/FTY_2T/Data/20250920_Xu_zfpl1l2_vs_l1")
getwd()

library(Seurat)
options(future.globals.maxSize = 1e10)

#load data
#读取10x数据，data.dir参数指定存放文件的路径

zfpl1l2_data <- Read10X(data.dir = "/media/dell/FTY_2T/Data/20250717_Xu_E13_ZFP36l1_l2_dcko_analysis/source_data/E20250449-01-01_RNA_中间文件/Emxl-L1L2/filtered_feature_bc_matrix")

#创建Seurat对象

zfpl1l2 <- CreateSeuratObject(counts = zfpl1l2_data,
                                 project = "zfpl1l2_data",
                                 min.features = 200,
                                 min.cells = 3)

#读取10x数据，data.dir参数指定存放文件的路径

zfpl1_data <- Read10X(data.dir = "/media/dell/FTY_2T/Data/20250717_Xu_E13_ZFP36l1_l2_dcko_analysis/source_data/E20250449-01-01_RNA_中间文件/20241206_Zfp36l1/E20240712-bc-01-02(RNA matrix)/E20240712-bc-01-02(RNA matrix)/zfp3611")

#创建Seurat对象

zfpl1 <- CreateSeuratObject(counts = zfpl1_data,
                              project = "zfpl1_data",
                              min.features = 200,
                              min.cells = 3)

zfpl1l2$Sample<-"l1l2"
zfpl1$Sample<-"l1"

###########QC############
data.combined<-merge(zfpl1,zfpl1l2)
table(data.combined$Sample)
grep(pattern = "^mt-",rownames(data.combined),value = TRUE) #返回具体基因名（小鼠）

data.combined[["percent.mt"]] <- PercentageFeatureSet(data.combined, pattern = "^mt-")
VlnPlot(data.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)

data.combined <- subset(data.combined, subset = Sample =="l1l2" & nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 2 |
                          Sample =="l1" & nFeature_RNA > 1200 & nFeature_RNA < 4500 & percent.mt < 2)
table(data.combined$Sample)

##########整合方法的选择############

data.combined <- NormalizeData(data.combined)
data.combined <- FindVariableFeatures(data.combined)
data.combined <- ScaleData(data.combined)
data.combined <- RunPCA(data.combined)

data.combined <- IntegrateLayers(
  object = data.combined, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

data.combined <- IntegrateLayers(
  object = data.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

data.combined[["RNA"]] <- JoinLayers(data.combined[["RNA"]])

##########可视化############

data.combined <- RunUMAP(data.combined, reduction = "integrated.cca", dims = 1:30,reduction.name = "umap_cca",
                         min.dist = 0.5,seed.use = 42)

#聚类
data.combined <- FindNeighbors(data.combined, reduction = "integrated.cca",dims = 1:30)
data.combined <- FindClusters(data.combined, random.seed = 0, resolution = 0.5,
                              cluster.name = "cca")

p1<-DimPlot(data.combined,split.by = "Sample",group.by = "cca",
            reduction = "umap_cca",label=T)

#########可视化
data.combined <- RunUMAP(data.combined, reduction = "integrated.rpca", dims = 1:30,reduction.name = "umap_rpca",
                         min.dist = 0.5,seed.use = 42)

#聚类
data.combined <- FindNeighbors(data.combined, reduction = "integrated.rpca",dims = 1:30)
data.combined <- FindClusters(data.combined, random.seed = 0, resolution = 0.5,
                              cluster.name = "rpca")

p2<-DimPlot(data.combined,split.by = "Sample",group.by = "rpca",
            reduction = "umap_rpca",label=T)


p1 / p2
