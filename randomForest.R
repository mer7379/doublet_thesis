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
library(fields)
#library(KernSmooth)
#library(parallel)
library(Seurat)
#library(BiocSingular)
# Load in the UMI matrix
library(randomForest)
pbmc.umis <- readRDS("pbmc_umi_mtx.rds")

### newww ##
#Load in the HTO count matrix
pbmc.htos <- readRDS("pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

#creat Seurat Object
sce.seu <- CreateSeuratObject(counts = pbmc.umis)

sce.hto <- CreateSeuratObject(counts = pbmc.htos)

#quality control preprocessing
sce.seu[["percent.mt"]] <- PercentageFeatureSet(sce.seu, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(sce.seu@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(sce.seu, features = c("nFeature_RNA",
                              "nCount_RNA", "percent.mt"), ncol = 3)

#Normalizing with log normalization
sce.seu <- NormalizeData(sce.seu, normalization.method = "LogNormalize", scale.factor = 10000)

sce.hto <- NormalizeData(sce.hto, normalization.method = "LogNormalize", scale.factor = 10000)

#Perform linear dimensional reduction -HTO
sce.hto <- HTODemux(sce.hto, assay = "RNA", positive.quantile = 0.99)
# First, we will remove negative cells from the object -HTO
#sce.seu <- subset(sce.seu, idents = "Negative", invert = TRUE)

#sce.seu <- RunPCA(sce.seu, features = VariableFeatures(object = sce.seu))
#sce.seu <- RunUMAP(sce.seu, dims = 1:10)

#find and scale variable features/ feature selection
sce.seu <- FindVariableFeatures(sce.seu, selection.method = "vst", nfeatures = 2000)
sce.seu <- ScaleData(sce.seu, features = VariableFeatures(sce.seu))

#library(solitude)
data.rf <- data.frame(sce.seu$nFeature_RNA,sce.seu$nCount_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$percent.mt)
data.rf$HTO <- as.factor(sce.hto$RNA_classification.global)

#random forest
set.seed(111)
rf.boston <- randomForest(data.rf$HTO~.,data=data.rf ,
                        mtry=4, importance =TRUE)
importance(rf.boston)
varImpPlot(rf.boston)

