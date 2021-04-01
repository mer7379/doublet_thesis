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
library(KernSmooth)
library(ROCR)
library(parallel)
library(DoubletFinder)
library(Seurat)
library(BiocSingular)
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
summary(sce.finder$pANN_0.25_0.01_1243)

##################################
######## visualisation ###########
##################################
plot_dl <- DimPlot(sce.finder,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = "DF.classifications_0.25_0.01_1243" )
plot_hto <- DimPlot(sce.finder,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = "HTO_classification.global")
plot_dl
plot_hto

##tsne plots
pbmc.tsne <- RunTSNE(sce.finder, dims = 1:8, perplexity = 100)
DimPlot(pbmc.tsne)
TSNEPlot(pbmc.tsne)

# select ground truth and pANN vectors from DF
mydata <- data.frame(sce.finder$HTO_classification.global,sce.finder$pANN_0.25_0.01_1243)

#split 70% train and 30%test
dt <- sort(sample(nrow(mydata),nrow(mydata)*0.7))
train <- mydata[dt,]
test <- mydata[-dt,]
test$y <- factor(test$sce.finder.HTO_classification.global, levels = c("Doublet","Singlet"), labels = c(1,0))
train$y<- factor(train$sce.finder.HTO_classification.global, levels = c("Doublet","Singlet"), labels = c(1,0))
#fit logistic regration
model_reg <- glm(train$y~train$sce.finder.pANN_0.25_0.01_1243, data= train, family = binomial(link = "logit"))
pred<-predict.glm(model_reg, newdata = data.frame(test$y,test$sce.finder.pANN_0.25_0.01_1243), type = c("response"))
summary(model_reg)

##############################################
#########Performance analysis using ROC#######
##############################################

#BiocManager::install("pROC")
library(pROC)
g<-roc(train$y~pred)
#plot.roc(train$y~pred, 
#         percent = TRUE,
 #        partial.auc=c(100, 90),
  #       partial.auc.correct=TRUE, print.auc=TRUE,
#        #display pAUC value on the plot with following options:
#         print.auc.pattern = "Corrected pAUC (100-90%% SP):\n%.1f%%",
#         print.auc.col = "#1c61b6",
#         auc.polygon = TRUE, 
#         auc.polygon.col = "#1c61b6",       # show pAUC as a polygon
#         max.auc.polygon = TRUE, 
#         max.auc.polygon.col = "#1c61b622", # also show the 100% polygon
#         main = "Partial AUC (pAUC)")
plot(g)
abline(a=0, b= 1)


