
# quality control
#setwd("/Users/ashenafiduba/Downloads/thesis/data")
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

#Load in the HTO count matrix
pbmc.htos <- readRDS("pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and
#HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

#setup single cell object
# Setup Seurat object and SingleCellExperiment
#pbmc.seu <- CreateSeuratObject(counts = pbmc.umis)

#pbmc.sce <- as.SingleCellExperiment(pbmc.seu)
#pbmc.sce
sce <- SingleCellExperiment(assays = list(counts = pbmc.umis))
sce

###quality control

#QC 1 Retrieving the mitochondrial transcripts using genomic locations included in
# the row-level annotation for the SingleCellExperiment
location <- rowRanges(sce)
is.mito <- any(seqnames(location)=="MT")

# ALTERNATIVELY: using resources in AnnotationHub to retrieve chromosomal
# locations given the Ensembl IDs; this should yield the same result.
#library(AnnotationHub)
#library(ensembldb)
#ens.mm.v97 <- AnnotationHub()[["AH73905"]]
#chr.loc <- mapIds(ens.mm.v97, keys=rownames(sce),
 #                 keytype="GENEID", column="SEQNAME")
#is.mito.alt <- which(chr.loc=="MT")
#########                   #####################
# calculating the QC metrics using perCellQCMetrics()
library(scater)
df <- perCellQCMetrics(sce, subsets=list(Mito=which(is.mito)))
df
#or alternatively This computes and appends the per-cell QC statistics 
#to the colData of the SingleCellExperiment object, allowing us to 
#retain all relevant information in a single object for later manipulation.
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
colnames(colData(sce))

#identifying low quality cells with fixed threshholds
qc.lib <- df$sum < 1e5
qc.nexprs <- df$detected < 5e3
#qc.spike <- df$altexps_ERCC_percent > 10
qc.mito <- df$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.mito

DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito), Total=sum(discard))

#QC 2 With adaptive thresholds
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

attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")
##added from quality pdf

par(mfrow=c(1,2))
hist(qc.lib2/1e6, xlab="Library sizes (millions)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(qc.nexprs2, xlab="Number of expressed genes", main="",
     breaks=20, col="grey80", ylab="Number of cells")
#altertnatives
reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))

##checking diagnostic plots

sce$discard <- reasons$discard                                               
plotColData(sce, x="sum", y="subsets_Mito_percent", colour_by = "discard")
#### for the poster session
#removing low quality cells
sce1 <- sce[,!reasons$discard]

#--- normalization ---#
set.seed(11111)
library(scran)
#sce <- logNormCounts(sce1)
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

#--- clustering ---#
#library(tables)
library(SingleCellExperiment)
snn.gr <- buildSNNGraph(sce1, use.dimred="PCA", k=25)
colLabels(sce1) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

#doublet cell detection
library(BiocSingular)
set.seed(100)

# Setting up the parameters for consistency with denoisePCA();
# this can be changed depending on your feature selection scheme.
library(scran)
dbl.dens <- doubletCells(sce1, subset.row=top.mam,
                         d=ncol(reducedDim(sce1)))
summary(dbl.dens)
head(dbl.dens)
sce1$DoubletScore <- log10(dbl.dens+1)
plotTSNE(sce1, colour_by="DoubletScore")

plotColData(sce1, x="label", y="DoubletScore", colour_by="label")
