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
library(parallel)
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

library(solitude)

data1 <- data.frame(sce.seu$nFeature_RNA,sce.seu$nCount_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$percent.mt) #1021

data2 <- data.frame(sce.seu$nFeature_RNA,sce.seu$nCount_RNA) #806
data3 <- data.frame(sce.seu$nFeature_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA)#1116
data4 <- data.frame(sce.seu$nCount_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA
) #992
data5 <- data.frame(
  sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
  sce.seu$percent.mt) #1178

data6 <- data.frame(sce.seu$nCount_RNA, sce.seu$percent.mt) #1016

data7 <- data.frame(sce.seu$nCount_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$percent.mt) #998
data8 <- data.frame(sce.seu$nFeature_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$percent.mt) #1098
data9 <- data.frame(sce.seu$nFeature_RNA,
                   sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$nCount_RNA) #998


#dt_iso <- sort(sample(nrow(data1),nrow(data1)*0.7))
#train_iso <- data1[dt_iso,]
#test_iso <- data1[-dt_iso,]
##HTO
#dt_hto <- sort(sample(length(sce.hto$RNA_classification.global),length(sce.hto$RNA_classification.global)*0.7))
#train_hto <- h[dt_hto,]
#test_hto <- h[-dt_hto,]

set.seed(3333)
index <- sample(ceiling(nrow(data2) * 0.005))# 256 sample_size = length(index)
iso_new <- isolationForest$new(sample_size = 256, num_trees = 100, seed = 130)
#iso_130 <- isolationForest$new(sample_size = length(data), seed = 130)
# fit for attrition data 
#iso_new$fit(train)  # training data

iso_new$fit(data2)  # training data
#?isolationForest

score <- iso_new$scores
# Obtain anomaly scores
set.seed(111)
iso_new_pred <- iso_new$predict(data2) # test
iso_new_pred[order(anomaly_score, decreasing = TRUE)]
#anomaly_pred <- (iso$predict(data))
summary(iso_new_pred$anomaly_score)
############################
#anomaly_score_130 <- scores_data_130$anomaly_score
#standard_dev_130 <- sd(anomaly_score_130)
# visualization of anomalyscore with respect to ground truth
#boxplot(iso_new_pred$anomaly_score~ test_hto)
length(iso_new_pred$anomaly_score[iso_new_pred$anomaly_score>0.6])

#boxplot(iso_new_pred$anomaly_score   ~ sce.hto$RNA_classification.global)


# select ground truth and pANN vectors from DF
#mydata <- data.frame(sce.finder$HTO_classification.global,sce.finder$pANN_0.25_0.2_1243)

mydata <- data.frame(sce.hto$RNA_classification.global,iso_new_pred$anomaly_score)

#split 70% train and 30%test
set.seed(2222)
dt <- sort(sample(nrow(mydata),nrow(mydata)*0.7))
train <- mydata[dt,]
test <- mydata[-dt,]
test$y <- factor(test$sce.hto.RNA_classification.global, levels = c("Doublet","Singlet"), labels = c(1,0))
train$y<- factor(train$sce.hto.RNA_classification.global, levels = c("Doublet","Singlet"), labels = c(1,0))
#fit logistic regration
model_reg <- glm(train$y~train$iso_new_pred.anomaly_score, data= train, family = binomial(link = "logit"))
pred<-predict.glm(model_reg, newdata = data.frame(test$y,test$iso_new_pred.anomaly_score), type = c("response"))
summary(model_reg)

##############################################
#########Performance analysis using ROC#######
##############################################

#BiocManager::install("pROC")

library(pROC)
#g1<-roc(train$y~pred, percent=TRUE)
#g2<-roc(train$y~pred, percent=TRUE)

#roc1 <- plot.roc(train$y~pred, main="ROC comparison", percent=TRUE, col= "red", print.auc = TRUE)
#roc5 <- lines.roc(train$y~pred, percent=TRUE,  sp = seq(0, 100, 5), reuse.auc = TRUE, auc.polygon = TRUE)

ROC6 <- plot(roc(train$y, pred), print.auc = TRUE)#BLACK AUC 0.736
ROC9 <- plot(roc(train$y, pred), print.auc = TRUE, col ="green", print.auc.y = 0.4, add = TRUE)
ROC1 <- plot(roc(train$y, pred), print.auc = TRUE, col ="RED", print.auc.y = 0.3, add = TRUE)
ROC2 <- plot(roc(train$y, pred), print.auc = TRUE, col ="ORANGE", print.auc.y = 0.35, add = TRUE)#AUC .772
ROC3 <- plot(roc(train$y, pred), print.auc = TRUE, col ="BLUE", print.auc.y = 0.2, add = TRUE)#AUC .772
ROC8 <- plot(roc(train$y, pred), print.auc = TRUE, col ="DARK GREEN", print.auc.y = 0.25, add = TRUE)#AUC .708
ROC5 <- plot(roc(train$y, pred), print.auc = TRUE, col ="DARK RED", print.auc.y = 0.45, add = TRUE)#AUC .551
ROC4 <- plot(roc(train$y, pred), print.auc = TRUE, col ="PINK", print.auc.y = 0.5, add = TRUE)#AUC .716 
ROC7 <- plot(roc(train$y, pred), print.auc = TRUE, col ="PURPLE", print.auc.y = 0.55, add = TRUE)#AUC 0.701 
 # the highest AUC is obtained by using data 2

#plot(g2, add = TRUE)
#g
#plot.roc(g1, g2)
#plot(g1, g2)
#abline(a=0, b= 1)
#boxplot(mydata$iso_new_pred.anomaly_score~ mydata$sce.hto.RNA_classification.global)


#### confusion matrix
library(caret)
library(e1071)
m <- confusionMatrix(factor(mydata$iso_new_pred.anomaly_score)
                     ,factor(mydata$sce.hto.RNA_classification.global))
m

