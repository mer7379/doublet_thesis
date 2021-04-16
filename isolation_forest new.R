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

library(solitude)
data <- data.frame(sce.seu$nFeature_RNA,sce.seu$nCount_RNA,
               sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
              sce.seu$percent.mt)

set.seed(1)
dt <- sample(ceiling(nrow(data) * 0.5))
train <- data[dt,]
test <- data[-dt,]


# initiate an isolation forest
iso <- isolationForest$new(sample_size = length(train))
# fit for attrition data 
iso$fit(train)  # training data
# Obtain anomaly scores
scores_train <- iso$predict(train)
scores_train[order(anomaly_score, decreasing = TRUE)]
# predict scores for test data (50% sample)
scores_unseen <- iso$predict(test)
scores_unseen[order(anomaly_score, decreasing = TRUE)]
#fit.train <- as.factor(ifelse(scores_train$anomaly_score >=0.70, 1, 0)) #1=doublet, 0=singlet
#pred.test <- as.factor(ifelse(scores_unseen$anomaly_score >=0.70, 1, 0))

#######################################
################performance############
          ################ 
set.seed(111)
library(pROC) # auc=0.6818
set.seed(111)
g<-roc(scores_train$anomaly_score~scores_unseen$anomaly_score) 
plot(g)

abline(a=0, b= 1)
g

#anomaly detection
# quantiles of anomaly scores
quantile(scores_unseen$anomaly_score
         , probs = seq(0.5, 1, length.out = 11))

################## by taking two variables only ##########
##########         : nFeature/nCount and percent.mt

#create isolation forest using isolationForest function from solitude package with default parameters

iforest<- isolationForest$new()
data <- data.frame(sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$percent.mt)
set.seed(1)
dt <- sort(sample(nrow(data),nrow(data)*0.5))
train <- data[dt,]
test <- data[-dt,]

iforest$fit(train)
# Obtain anomaly scores
scores_train <- iforest$predict(train)
scores_train[order(anomaly_score, decreasing = TRUE)]

#predict outliers within dataset
scores_unseen <- iforest$predict(test)
score_outlier <- as.factor(ifelse(scores_unseen$anomaly_score >= 0.70, "doublet", "singlet")) #i choose >0.7 as doublet
score_outlier_train <- as.factor(ifelse(scores_train$anomaly_score >= 0.70, "doublet", "singlet"))

#plot data again

Var1 <- test$sce.seu.nFeature_RNA.sce.seu.nCount_RNA
Var2 <- test$sce.seu.percent.mt
library(ggplot2)
ggplot(test, aes(x = Var1, y = Var2, color = score_outlier)) + 
  geom_point(shape = 1, alpha = 0.5) +
  labs(x = "x", y = "y") +
  labs(alpha = "", colour="Legend")

library(pROC) # result ==> auc=525
g<-roc(scores_train$anomaly_score~scores_unseen$anomaly_score) 
#g<-roc(scores_train$anomaly_score~scores_unseen$anomaly_score) 

plot(g)
abline(a=0, b= 1)
g


################## by taking two variables only ##########
##########         : nFeature/nCount and nCount_RNA
#create isolation forest using isolationForest function from solitude package with default parameters

iforest<- isolationForest$new()
data <- data.frame(sce.seu$nFeature_RNA/sce.seu$nCount_RNA,
                   sce.seu$nCount_RNA)
set.seed(1)
dt <- sort(sample(nrow(data),nrow(data)*0.5))
train <- data[dt,]
test <- data[-dt,]

iforest$fit(train)
# Obtain anomaly scores for the train dataset
scores_train <- iforest$predict(train)
scores_train[order(anomaly_score, decreasing = TRUE)]

#predict outliers within dataset
scores_unseen <- iforest$predict(test)
score_outlier <- as.factor(ifelse(scores_unseen$anomaly_score >= 0.70, "doublet", "singlet")) #i choose >0.7 as doublet
score_outlier_train <- as.factor(ifelse(scores_train$anomaly_score >= 0.70, "doublet", "singlet"))

#plot data 

Var1 <- test$sce.seu.nFeature_RNA.sce.seu.nCount_RNA
Var2 <- test$sce.seu.nCount_RNA
library(ggplot2)
ggplot(test, aes(x = Var1, y = Var2, color = score_outlier)) + 
  geom_point(shape = 1, alpha = 0.5) +
  labs(x = "x", y = "y") +
  labs(alpha = "", colour="Legend")

library(pROC) # result ==> auc=0.5147
set.seed(3333)
g<-roc(scores_unseen$anomaly_score~scores_train$anomaly_score) 
plot(g)
abline(a=0, b= 1)
g
