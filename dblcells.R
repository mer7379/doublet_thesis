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

library(scater)
library(SingleCellExperiment)
library(Seurat)

# Load in the UMI matrix
pbmc.umis <- readRDS("pbmc_umi_mtx.rds")
sce <- SingleCellExperiment(assays = list(counts = pbmc.umis))
sce

#### QUALITY CONTROL (QC) ###

# Retrieving the mitochondrial transcripts using genomic locations included in
# the row-level annotation for the SingleCellExperiment
location <- rowRanges(sce)
is.mito <- any(seqnames(location)=="MT")

#########                   #####################
# calculating the QC metrics using perCellQCMetrics()
df <- perCellQCMetrics(sce, subsets=list(Mito=which(is.mito)))
df
#or alternatively This computes and appends the per-cell QC statistics 
#to the colData of the SingleCellExperiment object, allowing us to 
#retain all relevant information in a single object for later manipulation.
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
colnames(colData(sce))

#QC 1 identifying low quality cells with fixed threshholds
qc.lib <- df$sum < 1e5
qc.nexprs <- df$detected < 5e3
#qc.spike <- df$altexps_ERCC_percent > 10 
qc.mito <- df$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.mito

DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito), Total=sum(discard))

#QC 2 With adaptive thresholds: remove cells that are outliers
# / this method is selected
qc.lib2 <- isOutlier(df$sum, log=TRUE, type="lower")
qc.nexprs2 <- isOutlier(df$detected, log=TRUE, type="lower")

attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")

qc.mito2 <- isOutlier(df$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")

discard2 <- qc.lib2 | qc.nexprs2 | qc.mito2
table(discard2)
# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2),
          MitoProp=sum(qc.mito2), Total=sum(discard2))

#attr(qc.lib2, "thresholds")
#attr(qc.nexprs2, "thresholds")

#altertnatives
reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))

#removing low quality cells
sce1 <- sce[,!reasons$discard]

#--- normalization ---#
set.seed(11111)
library(scran)
sce1 <- logNormCounts(sce1)

#--- variance-modelling ---#
set.seed(11111)
dec.mam <- modelGeneVarByPoisson(sce1)
top.mam <- getTopHVGs(dec.mam, prop=0.1)

#--- dimensionality-reduction ---#
library(BiocSingular)
set.seed(11111)
sce1 <- denoisePCA(sce1, technical=dec.mam, subset.row=top.mam)
sce1 <- runTSNE(sce1, dimred="PCA")
ncol(reducedDim(sce1))

#--- clustering ---#
#library(tables)
set.seed(1000)
library(SingleCellExperiment)
snn.gr <- buildSNNGraph(sce1, use.dimred="PCA", k=25)
colLabels(sce1) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(colLabels(sce1))

#___ visualisation ___#
set.seed(1000)
sce1 <- runTSNE(sce1, use_dimred="PCA")
plotTSNE(sce1, colour_by="label")

#___   doublet cell detection ___#
set.seed(1111)
library(scran)
dbl.dens <- doubletCells(sce1, subset.row=top.mam, 
                         d=ncol(reducedDim(sce1)))
summary(dbl.dens)
sce1$DoubletScore <- log10(dbl.dens+1)
plotTSNE(sce1, colour_by="DoubletScore")

### similar result with the above 
sce1$DoubletScore1 <- dbl.dens
plotTSNE(sce1, colour_by="DoubletScore")
plotColData(sce1, x="label", y="DoubletScore", colour_by="label")


