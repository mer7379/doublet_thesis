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

#quality control preprocessing
sce.seu[["percent.mt"]] <- PercentageFeatureSet(sce.seu, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(sce.seu@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(sce.seu, features = c("nFeature_RNA",
                              "nCount_RNA", "percent.mt"), ncol = 3)

#Normalizing with log normalization
sce.seu <- NormalizeData(sce.seu, normalization.method = "LogNormalize", scale.factor = 10000)

#find and scale variable features/ feature selection
sce.seu <- FindVariableFeatures(sce.seu, selection.method = "vst", nfeatures = 2000)
sce.seu <- ScaleData(sce.seu, features = VariableFeatures(sce.seu))

#Adding HTO data as an independent assay
# Add HTO data as a new assay independent from RNA
sce.seu[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
Assays(sce.seu)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
sce.seu <- NormalizeData(sce.seu, assay = "HTO", normalization.method = "CLR")

#####################################
###############################
#Perform linear dimensional reduction
sce.seu <- HTODemux(sce.seu, assay = "HTO", positive.quantile = 0.99)
# First, we will remove negative cells from the object
sce.seu <- subset(sce.seu, idents = "Negative", invert = TRUE)

sce.seu <- RunPCA(sce.seu, features = VariableFeatures(object = sce.seu))
sce.seu <- RunUMAP(sce.seu, dims = 1:10)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(sce.seu, PCs = 1:10, sct = FALSE)
gt.calls <- sce.seu@meta.data[rownames(sweep.res.list[[1]]), "HTO_classification.global"]  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats<- summarizeSweep(sweep.res.list, GT = TRUE, GT.calls = gt.calls)
bcmvn <- find.pK(sweep.stats)

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

### visualisation
plot_dl <- DimPlot(sce.finder,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = "DF.classifications_0.25_0.01_1243" )
plot_hto <- DimPlot(sce.finder,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = "HTO_classification.global")


#pbmc.tsne <- RunTSNE(sce.finder, dims = 1:8, perplexity = 100)
#DimPlot(pbmc.tsne)
#next visualization
########## clustering
#sce.seu <- FindNeighbors(sce.seu, dims = 1:10)
#sce.seu <- FindClusters(sce.seu, resolution = 0.5)
#look at cluster IDs of the first 6 cells
#head(Idents(sce.seu), 6)

