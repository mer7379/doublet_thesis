setwd("/Users/ashenafiduba/Downloads/thesis/data")
#install.packages("BiocManager")
## the command below is a one-line shortcut for:
## library(BiocManager)
## install("SingleCellExperiment")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("Seurat")
#BiocManager::install(c('scater', 'scran', 'uwot'))
#BiocManager::install("AnnotationHub")
#BiocManager::install("ensembldb")
#BiocManager::install("DoubletFinder")
#BiocManager::install("doubletCells")
#BiocManager::install("DoubletDecon")
#library(scater)
#library(SingleCellExperiment)
library(Seurat)
# Load in the UMI matrix
pbmc.umis <- readRDS("pbmc_umi_mtx.rds")

#creat Seurat Object
sce.seu <- CreateSeuratObject(counts = pbmc.umis)
#quality control preprocessing
sce.seu[["percent.mt"]] <- PercentageFeatureSet(sce.seu, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(sce.seu@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(sce.seu, features = c("nFeature_RNA",
                        "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

#plot1 <- FeatureScatter(sce.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(sce.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
sce.seu <- subset(sce.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalizing with log normalization
sce.seu <- NormalizeData(sce.seu, normalization.method = "LogNormalize", scale.factor = 10000)

#find and scale variable features/ feature selection

sce.seu <- FindVariableFeatures(sce.seu, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce.seu), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce.seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#sce.seu <- FindVariableFeatures(sce.seu, selection.method = "mean.var.plot")
#sce.seu <- ScaleData(sce.seu, features =  VariableFeatures(sce.seu))

##Scaling the data
all.genes <- rownames(sce.seu)
sce.seu <- ScaleData(sce.seu, features = all.genes)

#Perform linear dimensional reduction
sce.seu <- RunPCA(sce.seu, features = VariableFeatures(object = sce.seu))
#visualisation
VizDimLoadings(sce.seu, dims = 1:2, reduction = "pca")
DimHeatmap(sce.seu, dims = 1, cells = 500, balanced = TRUE)

#clustering
sce.seu <- RunUMAP(sce.seu, dims = 1:10)

##### DoubletFinder Function ###########
#BiocManager::install("DoubletFinder")
library(fields)
library(KernSmooth)
library(ROCR)
library(parallel)
library(DoubletFinder)

nExp_poi <- round(0.075*nrow(sce.seu@meta.data))
sce.finder <- doubletFinder_v3(sce.seu, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
summary(sce.finder)

########## clustering
#sce.seu <- FindNeighbors(sce.seu, dims = 1:10)
#sce.seu <- FindClusters(sce.seu, resolution = 0.5)
#look at cluster IDs of the first 6 cells
#head(Idents(sce.seu), 6)




