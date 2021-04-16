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

sce.dbl <- as.SingleCellExperiment(sce.seu)
set.seed(555)
sce.dbl <- runTSNE(sce.dbl, dimred="PCA")
sce_bcds <- bcds(sce.dbl)
#sce_bcds$bcds_score
#- Doublet scores are now available via colData
CD_bcds  <- colData(sce_bcds)
#head(cbind(CD_cxds$cxds_score))
######## visualisation ###########
##################################
##tsne plots
sce.dbl$doubletScore <- CD_bcds$bcds_score
plotTSNE(sce.dbl, colour_by="doubletScore")

# select ground truth and pANN vectors f
mydata <- data.frame(sce_bcds$HTO_classification.global,sce_bcds$bcds_score)

#split 70% train and 30%test
dt <- sort(sample(nrow(mydata),nrow(mydata)*0.7))
train <- mydata[dt,]
test <- mydata[-dt,]
test$y <- factor(test$sce_bcds.HTO_classification.global, levels = c("Doublet","Singlet"), labels = c(1,0))
train$y<- factor(train$sce_bcds.HTO_classification.global, levels = c("Doublet","Singlet"), labels = c(1,0))
#fit logistic regration
model_reg <- glm(train$y~train$sce_bcds.bcds_score, data= train, family = binomial(link = "logit"))
pred<-predict.glm(model_reg, newdata = data.frame(test$y,test$sce_bcds.bcds_score), type = c("response"))
summary(model_reg)

##############################################
#########Performance analysis using ROC#######
##############################################

library(pROC)
set.seed(4444)
g<-roc(train$y~pred) #auc=0.8099
plot(g)
abline(a=0, b= 1)
g




