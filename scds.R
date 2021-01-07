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
#BiocManager::install("scds")
#library(scater)
library(SingleCellExperiment)
#library(Seurat)
library(scds)
# Load in the UMI matrix
pbmc.umis <- readRDS("pbmc_umi_mtx.rds")
sce <- SingleCellExperiment(assays = list(counts = pbmc.umis))
sce

#- Annotate doublet using co-expression based doublet scoring:
sce_cxds <- cxds(sce)
#- Doublet scores are now available via colData
CD_cxds  <- colData(sce_cxds)
head(cbind(CD_cxds$cxds_score))

#- Annotate doublet using binary classification based doublet scoring:
sce_bcds <- bcds(sce)
#rm(sce_bcds)
#- Doublet scores are now available via colData:
CD_bcds  <- colData(sce_bcds)
head(cbind(CD_bcds$bcds_score))
rm(CD_bcds)


